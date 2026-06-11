import json

def convert_id(id_str: str) -> str | None:
    if id_str is None:          # ← guard here, not in convert_node
        return None
    if '_' in id_str:
        parts = id_str.rsplit('_', 1)
        return f"{parts[0]}:{parts[1]}"
    return id_str

def convert_node(node: dict) -> dict | None:
    node_id = convert_id(node.get("id"))
    if node_id is None:
        return None
    return {
        "id": node_id,
        "label": node["label"],
        "metadata": node["metadata"],
        "children": [c for c in (convert_node(ch) for ch in node.get("children", [])) if c is not None]
    }

def convert_file(input_path: str, output_path: str):
    with open(input_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    converted = [c for c in (convert_node(node) for node in data) if c is not None]

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump({"data": converted}, f, indent=2)

    print(f"Done. Written to {output_path}")

if __name__ == "__main__":
    import sys
    input_path  = sys.argv[1] if len(sys.argv) > 1 else "input.json"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "output.json"
    convert_file(input_path, output_path)
