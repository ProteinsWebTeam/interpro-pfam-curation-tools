#!/usr/bin/env python3


import sys
import cx_Oracle
import traceback
import re


class pfam_duf:
    def __init__(self):
        self.connection = None
        self.cursor = None
        self.list_duf = dict()
        self.count_interpro_dom = dict()

    def getConnection(self, user, password, schema):
        """
		Set database connection
		"""

        connectString = "".join([user, "/", password, "@", schema])
        try:
            # subprocess.run(["source", "~oracle/ora112setup.sh"])
            self.connection = cx_Oracle.connect(connectString)
            self.cursor = self.connection.cursor()
        except:
            stackTrace = traceback.format_exc()
            print(stackTrace)
            if "invalid username" in stackTrace:
                print("Could not connect to {0} as user {1}".format(schema, user))
                print(
                    "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
                )
                print('eg "ops$craigm" or ops\$craigm')
                print("Otherwise the shell will remove the '$' and all subsequent characters!")
                sys.exit(1)

    def close_connection(self):
        """
		Close database connection
		"""

        self.cursor.close()
        self.connection.close()

    def chunks(self, l, n):
        """Yield chunks of size n from iterable.

		Args:
			l (Iterable): The iterable to chunk
			n (int): Maximum number of items in chunk

		Yields:
			Iterable: Tuples of length n or less (final bucket will have only the remaining items)

		"""
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i : i + n]

    def get_pfam_duf_list(self):
        request = "select m.method_ac, m.name, e2m.ENTRY_AC, e.NAME \
                    from interpro.method m \
                    left join interpro.entry2method e2m on e2m.method_ac=m.method_ac \
                    left join interpro.entry e on e2m.entry_ac=e.entry_ac \
                    where m.method_ac like 'PF%' and m.name like 'DUF%' "

        self.cursor.execute(request)

        for row in self.cursor:
            self.count_interpro_dom[row[2]] = 0
            self.list_duf[row[0]] = {"dufid": row[1], "ipr": row[2], "name": row[3]}

    def get_nb_domain_per_interpro(self):
        print("Searching for average number of domains per InterPro entry")
        ipr_chunks = list(self.chunks(list(self.count_interpro_dom.keys()), 1000))

        for chunk in ipr_chunks:
            list_ipr_quote = [f"'{ipr}'" for ipr in chunk]
            request = f"select e2m.entry_ac, round(avg(sub.cpt)) \
                    from entry2method e2m \
                    join match m on m.method_ac=e2m.method_ac \
                    join ( \
                    select m.protein_ac, count(distinct method_ac) as cpt \
                    from INTERPRO.MATCH m \
                    join INTERPRO.PROTEIN p on p.protein_ac=m.protein_ac \
                    where m.DBCODE='H' and p.DBCODE='S' and p.FRAGMENT='N' \
                    group by m.protein_ac \
                    ) sub on sub.protein_ac = m.protein_ac \
                    where e2m.entry_ac in ({','.join(list_ipr_quote)}) and m.DBCODE='H' \
                    group by e2m.ENTRY_AC"

            # print(request)
            self.cursor.execute(request)
            for row in self.cursor:
                self.count_interpro_dom[row[0]] = row[1]

    def save_duf_list_in_file(self, outfile, outfileipr):

        with open(outfile, "w") as outf, open(outfileipr, "w") as outfipr:
            outf.write("pfamid,dufid,ipr,avg nb domains per ipr,ipr entry name\n")
            outfipr.write("pfamid,dufid,ipr,avg nb domains per ipr,ipr entry name\n")
            for pfamid, content in self.list_duf.items():
                ipr = content["ipr"]
                ipr_count = self.count_interpro_dom[ipr]

                outf.write(f"{pfamid},{content['dufid']},{ipr},{ipr_count},{content['name']}\n")
                if (
                    content["name"] != None
                    and not re.search(r"DUF", content["name"])
                    and not re.search(r"Uncharacterised", content["name"])
                ):
                    outfipr.write(
                        f"{pfamid},{content['dufid']},{ipr},{ipr_count},{content['name']}\n"
                    )

