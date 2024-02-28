import networkx as nx
import matplotlib.pyplot as plt

# First, let's load the JSON file to understand its structure
import json

file_path = 'LDLR_clusters-all-MiniLM-L6-v2.json'

# Load the JSON content
with open(file_path, 'r') as file:
    data = json.load(file)


# Create a new graph
G = nx.Graph()

# Add nodes and edges based on the JSON data structure
for group in data:
    group_id = group['id']  # Group ID
    G.add_node(group_id, label=f"Group {group_id}", color='lightblue')  # Add group node
    
    for item in group['items']:
        item_id = item['id']  # Item ID
        item_content = item['content']  # Item content
        G.add_node(item_id, label=item_content, color='lightgreen')  # Add item node
        G.add_edge(group_id, item_id)  # Connect group node to item node

# Drawing the network with labels and custom colors
pos = nx.spring_layout(G)  # Positioning the nodes using the spring layout
colors = nx.get_node_attributes(G, 'color').values()
labels = nx.get_node_attributes(G, 'label')

nx.draw(G, pos, with_labels=False, node_color=list(colors), alpha=0.7, node_size=700)
nx.draw_networkx_labels(G, pos, labels, font_size=8)

plt.title("Network Visualization of Groups and Items")
plt.axis('off')  # Turn off the axis
plt.show()

