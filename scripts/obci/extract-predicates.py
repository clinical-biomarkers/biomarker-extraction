import re
import json
from collections import defaultdict

def extract_predicates(owl_file_path):
    with open(owl_file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Pass 1: find all property URIs actually used in restrictions
    used_props = set(re.findall(
        r'<owl:onProperty rdf:resource="([^"]+)"',
        content
    ))

    # Find all ObjectProperty blocks
    prop_pattern = re.compile(
        r'<owl:ObjectProperty[^>]*>.*?</owl:ObjectProperty>',
        re.DOTALL
    )

    ALLOWED_PREFIXES = ('OBCI:', 'RO:')
    result = {}

    for match in prop_pattern.finditer(content):
        block = match.group()

        xref = re.search(r'<oboInOwl:hasDbXref[^>]*>([^<]+)</oboInOwl:hasDbXref>', block)
        oid  = re.search(r'<oboInOwl:id[^>]*>([^<]+)</oboInOwl:id>', block)

        if xref and oid:
            # OBCI-style: both annotations present
            key = xref.group(1).strip()
            val = oid.group(1).strip()
        else:
            # RO-style: derive key from rdf:about URI, label as value
            about = re.search(r'<owl:ObjectProperty[^>]+rdf:about="([^"]+)"', block)
            label = re.search(r'<rdfs:label[^>]*>([^<]+)</rdfs:label>', block)
            if not about or not label:
                continue
            uri = about.group(1)
            ro_id = re.search(r'/(RO)_(\d+)$', uri)
            if not ro_id:
                continue
            # Restrict to RO properties actually used in OBCI axioms
            if uri not in used_props:
                continue
            key = f"{ro_id.group(1)}:{ro_id.group(2)}"
            val = label.group(1).strip()

        if key.startswith(ALLOWED_PREFIXES):
            result[key] = val

    return {"ontology": {"predicate": result}}

if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "ontology.owl"
    output = extract_predicates(path)
    print(json.dumps(output, indent=2))
