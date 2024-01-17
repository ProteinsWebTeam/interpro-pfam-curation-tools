import os
import pickle
import re
import sys
from argparse import ArgumentParser
from datetime import datetime

import cx_Oracle


PRONTO = "http://pronto.ebi.ac.uk:5000/search/?q="
LETTERS = (
    ('Alpha', 'Α', 'α'), 
    ('Beta', 'Β', 'β'), 
    ('Gamma', 'Γ', 'γ'), 
    ('Delta', 'Δ', 'δ'), 
    ('Epsilon', 'Ε', 'ε'), 
    ('Zeta', 'Ζ', 'ζ'), 
    ('Eta', 'Η', 'η'), 
    ('Theta', 'Θ', 'θ'), 
    ('Iota', 'Ι', 'ι'), 
    ('Kappa', 'Κ', 'κ'), 
    ('Lambda', 'Λ', 'λ'), 
    # ('Mu', 'Μ', 'μ'), 
    # ('Nu', 'Ν', 'ν'), 
    # ('Xi', 'Ξ', 'ξ'), 
    ('Omicron', 'Ο', 'ο'), 
    # ('Pi', 'Π', 'π'), 
    ('Rho', 'Ρ', 'ρ'), 
    ('Sigma', 'Σ', 'σ'), 
    ('Tau', 'Τ', 'τ'), 
    ('Upsilon', 'Υ', 'υ'), 
    ('Phi', 'Φ', 'φ'), 
    ('Chi', 'Χ', 'χ'), 
    ('Psi', 'Ψ', 'ψ'), 
    ('Omega', 'Ω', 'ω')
)

MAP_LETTERS = {
    ascii_letter.lower(): greek_lower
    for ascii_letter, greek_upper, greek_lower in LETTERS
}
GREEK_LETTERS = "".join(greek_lower for ascii_letter, greek_upper, greek_lower in LETTERS)

def ascii2greek(match):
    ascii_letter = match.group(1)
    greek_sub = MAP_LETTERS[ascii_letter.lower()]
    return match.group(0).replace(ascii_letter, greek_sub)

def replace_greek(parts_in):
    parts_out = []
    for item in parts_in:
        try:
            greek_sub = MAP_LETTERS[item.strip().lower()]
            item = item.replace(item.strip(), greek_sub)
        except KeyError:
            if item != "":
                print(f"Error greek letter not found: {item}")
        else:
            parts_out.append(item)
    return parts_out

class Finder:
    def __init__(self, to_replace):
        self.to_replace = to_replace

    def sub(self, match):
        # ascii_letter = match.group(1)
        entire_match = match.group(0)
        striped_entire_match = entire_match.strip()

        if striped_entire_match in self.to_replace:
            # print(f"before: '{entire_match}'")
            for letter in MAP_LETTERS:
                if re.search(letter, entire_match.lower()):
                    greek_letter = MAP_LETTERS[letter]
                    entire_match = entire_match.replace(letter, greek_letter)
            
            # print(f"after: '{entire_match}'")
            return entire_match
        else:
            return entire_match


def main():
    parser = ArgumentParser()
    parser.add_argument("--replace", help="file of matches to replace")
    parser.add_argument("--restore", help="file of abstracts to restore")
    parser.add_argument("--print", action="store_true", 
                        help="print new abstracts instead of updating the DB")
    args = parser.parse_args()

    to_replace = set()
    if args.replace:
        with open(args.replace, "rt") as fh:
            for line in fh:
                to_replace.add(line.rstrip())

    try:
        url = os.environ["INTERPRO_URL"]
    except KeyError:
        parser.error("INTERPRO_URL environment variable not set")
    
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    if args.restore:
        with open(args.restore, "rb") as fh:
            params = pickle.load(fh)
        
        s = input(f"Restore {len(params)} CABs [y/N]? ")
        if s not in ('y', 'Y'):
            cur.close()
            con.close()
            sys.stderr.write("Aborted\n")
            return
        
        cur.executemany(
            """
            UPDATE INTERPRO.COMMON_ANNOTATION
            SET TEXT = :1, COMMENTS = :2
            WHERE ANN_ID = :3
            """, params
        )

        con.commit()
        cur.close()
        con.close()
        sys.stderr.write("Done\n")
        return
    
    cur.execute("SELECT ANN_ID, COMMENTS, TEXT FROM INTERPRO.COMMON_ANNOTATION")
    now = datetime.now()
    new_comment = f"Annotation updated automatically on {now:%Y-%m-%d %H:%M:%S}"
    backup = []
    params = []
    cases = {}

    finder = Finder(to_replace)
    for ann_id, comment, text in cur:
        new_text = text
        
        for ascii_letter, greek_upper, greek_lower in LETTERS:
            # pattern = rf"[^a-z]({ascii_letter}-?[a-z0-9\.\-]{{2,}})"
            """
            Letter followed by:
                - optional hypen + number(s)
                - optional hyphen or slash + text
            """
            #add greek_lower list to regex
            # pattern = rf"\b({ascii_letter})(?:-?\d+|-[a-z]{{2,}})\b"
            # pattern = rf"[ {GREEK_LETTERS}/+-]+({ascii_letter})(?:-?\d+|[ /+-]+[{GREEK_LETTERS}a-z\d/-]{{2,}})\b"
            # pattern = rf"[ /+-]+({ascii_letter})(?:-?\d+|[ /+-]+[{GREEK_LETTERS}a-z\d/-]{{2,}})\b"
            pattern = rf"([ \(/+-]+{ascii_letter})([ \(\)/+-]+[{ascii_letter}{GREEK_LETTERS}a-z\d/-]+)([ \(\)/+-]+[{ascii_letter}{GREEK_LETTERS}a-z\d/-]+)\b"

            for match in re.finditer(pattern, text, flags=re.I):
                entire_match = match.group(0)

                try:
                    cases[ascii_letter][entire_match] = ann_id
                except KeyError:
                    cases[ascii_letter] = {entire_match: ann_id}

                new_text = re.sub(pattern, finder.sub, new_text, flags=re.I)

        if new_text != text:
            backup.append((
                text,
                comment,
                ann_id
            ))
            params.append((
                new_text,
                new_comment,
                ann_id
            ))

    if to_replace:
        if not params:
            cur.close()
            con.close()
            sys.stderr.write("No CABs to update\n")
        elif args.print:
            cur.close()
            con.close()

            for i, (new_text, comment, ann_id) in enumerate(params):
                if i:
                    print('-' * 40)

                print(ann_id)
                print(new_text)
        else:
            s = input(f"Update {len(params)} CABs [y/N]? ")
            if s not in ('y', 'Y'):
                cur.close()
                con.close()
                sys.stderr.write("Aborted\n")
                return
                
            for param in params:
                try:
                    cur.execute(
                        """
                        UPDATE INTERPRO.COMMON_ANNOTATION
                        SET TEXT = :1, COMMENTS = :2
                        WHERE ANN_ID = :3
                        """, param
                    )
                except cx_Oracle.IntegrityError:
                    print(param)
            
            # cur.executemany(
            #     """
            #     UPDATE INTERPRO.COMMON_ANNOTATION
            #     SET TEXT = :1, COMMENTS = :2
            #     WHERE ANN_ID = :3
            #     """, params
            # )

            con.commit()
            
            cur.close()
            con.close()

            with open("abstracts_3.pickle", "wb") as fh:
                pickle.dump(backup, fh)
            
            sys.stderr.write("Done\n")
    else:
        cur.close()
        con.close()
        for ascii_letter in sorted(cases):
            matches = cases[ascii_letter]
            for match in sorted(matches):
                if not re.search('subunit', match):
                    print(f"{ascii_letter}\t{match}\t{PRONTO}{matches[match]}")


if __name__ == '__main__':
    main()

