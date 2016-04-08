from xml.dom.minidom import parse
import json

file = open("iJO1366.xml")
xml_dom = parse(file)

species = xml_dom.getElementsByTagName("species")
reactions = xml_dom.getElementsByTagName("reaction")

node_list = species + reactions

d = {}
for node in node_list:
    model_id = node.getAttribute("id")
    resources = node.getElementsByTagName("rdf:li")
    for r in resources:
        path = r.getAttribute("rdf:resource")
        if "kegg.compound" in path or "kegg.reaction" in path:
            kegg_id = path.split("/")[-1]
            d[model_id] = kegg_id

print(json.dumps(d, indent=4, sort_keys=True))