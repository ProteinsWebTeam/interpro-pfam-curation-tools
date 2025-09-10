-- connect to pronto database
-- select all DUF Pfam accessions from current release, add quotes around the Pfam accessions and separate them with commas. 
-- Add the Pfam accessions to the IN clause of the SQL queries below.

SELECT sp.signature_acc, p.accession, p.identifier, n.text
            FROM interpro.signature2protein sp
            INNER JOIN interpro.protein p
                ON sp.protein_acc = p.accession
            INNER JOIN interpro.protein2name pn
                ON pn.protein_acc = sp.protein_acc
            INNER JOIN interpro.protein_name n
                ON pn.name_id = n.name_id
            WHERE p.is_reviewed='t'
            AND lower(n.text) not like 'putative%'
            AND lower(n.text) not like '%domain-containing%'
            AND lower(n.text) not like '%uncharacterized%'
            AND lower(n.text) not like '%unknown%'
            AND sp.signature_acc IN ()
            order by sp.signature_acc;

SELECT p.accession, p.identifier, p.is_reviewed,
                   sp.signature_acc, ss.structure_id
            FROM interpro.signature2protein sp
            INNER JOIN interpro.signature2structure ss
                ON sp.signature_acc = ss.signature_acc
               AND sp.protein_acc = ss.protein_acc
            INNER JOIN interpro.protein p
                ON sp.protein_acc = p.accession
            WHERE sp.signature_acc IN ();
