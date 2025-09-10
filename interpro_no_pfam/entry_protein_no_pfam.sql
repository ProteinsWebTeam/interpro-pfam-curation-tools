
/* search for entries with proteins that are not matched by any pfam*/
SET SQLFORMAT csv

SPOOL entry_no_pfam.csv

    select distinct(e2m.entry_ac)
    from interpro.match m
    join interpro.entry2method e2m on m.method_ac=e2m.method_ac
    where dbcode!='H'
    and m.protein_ac not in 
    (
    select distinct(protein_ac)
    from interpro.match
    where dbcode='H'
    );

SPOOL OFF