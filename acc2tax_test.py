import acc2tax

acc2tax.init(16, "/mnt/data/nt/nucl_gb.accession2taxid.gz", 
	"/mnt/data/nt/taxdmp/nodes.dmp", 
	"/mnt/data/nt/taxdmp/names.dmp");
print(acc2tax.get_taxon("X52700.1"))
tax_id = acc2tax.get_taxon("X52700.1")

print(acc2tax.get_complete_taxonomy(tax_id))
print()
print(acc2tax.get_complete_taxonomy_dict(tax_id))
print()
print(acc2tax.get_complete_taxonomy_names_dict(tax_id))

tax_id = acc2tax.get_taxon("XR_002086441.1")
print(acc2tax.get_taxon("XR_002086441.1"))
print(acc2tax.get_complete_taxonomy_names_dict(tax_id))

acc = "X16634.1"
tax_id = acc2tax.get_taxon(acc)
print(acc2tax.get_taxon(acc))
print(acc2tax.get_complete_taxonomy_names_dict(tax_id))

tax_id = 5145
print(acc2tax.get_complete_taxonomy_names_dict(tax_id))
