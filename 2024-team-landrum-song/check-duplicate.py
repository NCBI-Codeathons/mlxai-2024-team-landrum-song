import json, sys
# https://chat.openai.com/share/1cef3a28-2555-4a79-9dc5-145cc2b1dcee 

gene = sys.argv[1]

# Load LDLR_test.json
with open('data/' + gene + '_variants_extracted.json', 'r') as file:
    ldlr_test_data = json.load(file)

# Load LDLR_clusters.json
with open(gene + '_clusters.json', 'r') as file:
    ldlr_clusters_data = json.load(file)

# Check the first few items of each to understand their structure
ldlr_test_data[:2], ldlr_clusters_data[:2]

# Step 2: Organize LDLR_clusters_data for efficient searching
# Exclude clusters with id: "-1" and create a mapping of content to cluster id

content_to_cluster_id = {}
for cluster in ldlr_clusters_data:
    if cluster['id'] != "-1":
        for item in cluster['items']:
            content_to_cluster_id[item['content'].lower()] = cluster['id']

# Step 3: Process LDLR_test_data to find records meeting the criteria
# For each variant_id, we will check if 2 or more text fields can be found in any content fields of clusters

# Organize LDLR_test_data by variant_id and collect texts
variant_id_to_texts = {}
for record in ldlr_test_data:
    variant_id = record['variant_id']
    text = record['condition']['text'].lower()
    if variant_id not in variant_id_to_texts:
        variant_id_to_texts[variant_id] = []
    variant_id_to_texts[variant_id].append(text)

# Find qualifying records
qualifying_records = []
total_reduced_records = 0 
counter_unique_variant_id_duplicate = 0  # count the number of unique variant_id with similar disease names 
for variant_id, texts in variant_id_to_texts.items():
    # Count texts found in any content fields of clusters
    found_texts_count = sum(1 for text in texts if text in content_to_cluster_id)
   
    # If 2 or more texts found, add all records with this variant_id to qualifying_records
    if found_texts_count >= 2:
        qualifying_records.extend([record for record in ldlr_test_data if record['variant_id'] == variant_id])
        total_reduced_records += found_texts_count -1 
        counter_unique_variant_id_duplicate = counter_unique_variant_id_duplicate + 1 

print('\nNumber of RCV records that can be potentially reduced: ', total_reduced_records)
print('Total nmber of RCV records: ', len(ldlr_test_data)) 

print('\nNumber of variants with groupable disease names: ', counter_unique_variant_id_duplicate)
print('Number of unique variants: ' , len(variant_id_to_texts)) # number of unique variant_id 
print()

# Let's check the number of qualifying records found
# len(qualifying_records), qualifying_records[:2]  # Display the count and first two for review

# Step 4: Save the qualifying records into a new JSON file

output_file_path = gene + '_qualifying_records.json'
with open(output_file_path, 'w') as outfile:
    json.dump(qualifying_records, outfile, indent=4)

print('Write ', output_file_path)



# # Adjusting the process to include only the qualifying matching records
# qualifying_matching_records = []
# for record in ldlr_test_data:
#     text = record['condition']['text'].lower()
#     # Check if the record's text is found in the content-to-cluster mapping
#     if text in content_to_cluster_id:
#         # Directly add this matching record to the list
#         qualifying_matching_records.append(record)
# 
# # Save the qualifying matching records into a new JSON file with indent=4
# output_file_path_matching = 'qualifying_matching_records.json'
# with open(output_file_path_matching, 'w') as outfile:
#     json.dump(qualifying_matching_records, outfile, indent=4)
# 
# 
# print('Write ', output_file_path_matching)

# # This time, we will check each record individually instead of adding all records with a matching variant_id
# individual_qualifying_records = []
# for record in ldlr_test_data:
#     variant_id = record['variant_id']
#     text = record['condition']['text'].lower()
#     # Check if the text of the current record is in any content fields of clusters
#     if text in content_to_cluster_id:
#         # Add this record to individual_qualifying_records
#         individual_qualifying_records.append(record)
# 
# # Save the individual qualifying records into a new JSON file with indent=4
# output_file_path_individual = 'individual_qualifying_records.json'
# with open(output_file_path_individual, 'w') as outfile:
#     json.dump(individual_qualifying_records, outfile, indent=4)


# print('Write ', output_file_path_individual)
