import sys
import networkx as nx
import matplotlib.pyplot as plt

# First, let's load the JSON file to understand its structure
import json

file_path = sys.argv[1]
#file_path = 'LDLR_clusters.json'
#file_path = 'SCN5A_clusters.json'
#file_path = 'KCNQ1_clusters.json'
#file_path = 'USH2A_clusters.json'
#file_path = 'TSC1_clusters.json'

# Load the JSON content
with open(file_path, 'r') as file:
    data = json.load(file)
# print(data)

from pyvis.network import Network

# Create a new interactive network
net = Network(notebook=True, height="1000px", width="100%")

# print(group[0])

n=0
group_id = 'A'

# Add nodes and edges based on the JSON data structure
for group in data:
    if group['id'] == '-1':
        group_id = 'Noise'
    else:
        group_id = chr(int(group['id']) +  65)  
    #print (group['id'])
    
    net.add_node(group_id, label=f"Group {group_id}", color='lightblue', size=40)
    
    for item in group['items']:
        # item_id = item['id']
        item_id = n 
        item_content = item['content']
        net.add_node(item_id, label=item_content, color='lightgreen', size=30)
        net.add_edge(group_id, item_id)
        n = n+1
    
    #group_id = chr(ord(group_id) + 1)

# Customize the network appearance
net.set_options("""
var options = {
  "nodes": {
    "font": {
      "size": 36
    }
  },
  "edges": {
    "color": {
      "inherit": true
    },
    "smooth": false
  },
  "physics": {
    "forceAtlas2Based": {
      "gravitationalConstant": -250,
      "centralGravity": 0.005,
      "springLength": 230,
      "springConstant": 0.18
    },
    "maxVelocity": 146,
    "solver": "forceAtlas2Based",
    "timestep": 0.35,
    "stabilization": {"iterations": 150}
  }
}
""")

# Display the network
net.show("interactive_" + file_path[0:-5] + ".html")


