import json

def check_json(file_path):
    with open(file_path, 'r') as file:
        raw_data = file.read()

    try:
        parsed_data = json.loads(raw_data)
        print("JSON validated successfully!")
    except json.JSONDecodeError as e:
        print(f"JSON validation failed: {e}")

# Call the function with your file
check_json("all-compare.json")
