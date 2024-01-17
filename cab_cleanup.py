#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script cleanups annotation blocks with white space and missing <p> blocks

@arguments [CONFIG_FILE]: file containing database connection and files location (see config.ini)

"""

import argparse
import os
import re
from configparser import ConfigParser
import traceback
import cx_Oracle
from datetime import datetime

def revert_annotation(user, password, schema):
    connectString = "".join([user, "/", password, "@", schema])

    try:
        # subprocess.run(["source", "~oracle/ora112setup.sh"])
        connection = cx_Oracle.connect(connectString)
        cur = connection.cursor()
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)

    date_old = datetime(2023, 2, 28, 23, 59, 00)
    date_new = datetime(2023, 3, 1, 00, 00, 00)
    date_last = datetime(2023, 3, 17, 00, 00, 00)
    request_test = """
        SELECT ann_id,timestamp
        FROM common_annotation_audit
        WHERE DBUSER!='INTERPRO' and timestamp < :1
        AND (ann_id, timestamp) in (
                SELECT ann_id, max(timestamp) 
                FROM common_annotation_audit 
                WHERE timestamp<:1
                GROUP BY ann_id
            )
    """
    request = """
        SELECT caa_old.ann_id, caa_old.len_old, caa_new.len_new, caa_old.text as text_old, caa_new.text as text_new, caa_new.timestamp
        FROM (
            SELECT ann_id, length(text) as len_old, text
            FROM common_annotation_audit 
            WHERE timestamp< :1 and DBUSER!='INTERPRO'
            AND (ann_id, timestamp) in (
                SELECT ann_id, max(timestamp) 
                FROM common_annotation_audit 
                WHERE timestamp< :2
                GROUP BY ann_id
            )
         ) caa_old
         JOIN (
             SELECT ann_id, length(text) as len_new, text, timestamp 
             FROM common_annotation_audit 
             WHERE timestamp> :3 AND timestamp< :4 AND DBUSER='INTERPRO'
         ) caa_new ON caa_old.ann_id=caa_new.ann_id
         WHERE caa_new.len_new < caa_old.len_old/2
     """

    update_cab = set()
    # cur.execute(request)
    for row in cur.execute(request, (date_old, date_old, date_new, date_last)):
        # print(row)
        ann_id = row[0]
        new_text = row[3]
        now = datetime.now()
        new_comment = f"Reverted automatic update on {now:%Y-%m-%d %H:%M:%S}"
        update_cab.add((new_text, new_comment, ann_id))

    print(len(update_cab))

    for item in update_cab:
        try:
            cur.execute("UPDATE COMMON_ANNOTATION SET TEXT=:1, COMMENTS=:2 WHERE ANN_ID=:3", item)
        except cx_Oracle.IntegrityError:
            print(item[2])
    
    connection.commit()

    cur.close()
    connection.close()
    
def delete_white_space(text) -> str:
    pattern_beg = r"\> +"
    pattern_end = r" +\<\/"

    skip = ["<reaction>", "<sub>", "<sup>", "<i>", "->"]
    found = 0
    for item in skip:
        if re.search(item, text):
            found += 1
    if found == 0:
        text = re.sub(pattern_beg, ">", text)
        text = re.sub(pattern_end, "</", text)
    
    return text

def clean_annotation(user, password, schema):
    connectString = "".join([user, "/", password, "@", schema])
    list_cdd = dict()

    try:
        # subprocess.run(["source", "~oracle/ora112setup.sh"])
        connection = cx_Oracle.connect(connectString)
        cursor = connection.cursor()
    except:
        stackTrace = traceback.format_exc()
        print(stackTrace)

    request = """
        SELECT ANN_ID, TEXT 
        FROM INTERPRO.COMMON_ANNOTATION
        """
    
    prog = re.compile(r"[a-z]\s{2,}[a-z]") #done check for space at the beginning of a cab
    prog2 = re.compile(r"^<p>[a-z]") #done add <p></p> tags
    prog3 = re.compile(r"([\w\W]+\w+)(\s{2,})(.+)") #search for more than 1 white space 
    prog4 = re.compile(r"<p>-") #done
    prog5 = re.compile(r"<p>(.+<ul>.+)</p>") #done
    prog6 = re.compile(r"^([\w\W]+)(\. </p>)([\n\w\W]+)") #done
    prog7 = re.compile(r" </li>\n\n")
    prog8 = re.compile(r"<p>\s+[a-zA-Z]")
    # prog9 = re.compile(r"\.\s+</li>")

    
    update_cab = []
    cursor.execute(request)
    _RE_COMBINE_WHITESPACE = re.compile(r"\s+")
    save_ann_id = ""
    count = 0
    for row in cursor:
        ann_id, text = row
        updated_text = delete_white_space(text)
        if updated_text != text:
            count += 1
            if count < 50:
                print(f"\n{ann_id}")
                print(text)
                print("\nnew:")
                print(updated_text)
                replace = input("Replace CAB?")
                if replace == "y":
                    now = datetime.now()
                    new_comment = f"Annotation updated automatically on {now:%Y-%m-%d %H:%M:%S}"
                    update_cab.append((updated_text, new_comment, ann_id))
            else:
                break
    #     for match in prog.findall(text):
    #         if save_ann_id != ann_id:
    #             count += 1
    #             if count < 100:
    #                 print(f"\n{ann_id}")
    #                 print(text)
    #                 updated_text = _RE_COMBINE_WHITESPACE.sub(" ", text).strip()
    #                 updated_text = updated_text.replace("</p> <p>", "</p>\n\n<p>")
    #                 print("\nnew:")
    #                 print(updated_text)
    #                 replace = input("Replace CAB?")
    #                 if replace == "y":
    #                     now = datetime.now()
    #                     new_comment = f"Annotation updated automatically on {now:%Y-%m-%d %H:%M:%S}"
    #                     update_cab.append((updated_text, new_comment, ann_id))
    #                 save_ann_id = ann_id
    #             else:
    #                 break
    # print(len(update_cab))
    
    for item in update_cab:
        try:
            cursor.execute("UPDATE COMMON_ANNOTATION SET TEXT=:1, COMMENTS=:2 WHERE ANN_ID=:3", item)
        except cx_Oracle.IntegrityError:
            print(item[2])
    
    connection.commit()

    cursor.close()
    connection.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    # database connection values
    user = config["database"]["user"]
    password = config["database"]["password"]
    schema = config["database"]["schema"]

    clean_annotation(user, password, schema)
    # revert_annotation(user, password, schema)