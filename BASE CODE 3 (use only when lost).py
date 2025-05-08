import streamlit as st
import requests
import networkx as nx
import plotly.graph_objects as go
import io

# -------------- Background Styling --------------
def style_app_container():
    css = f"""
    <style>
    .block-container {{
        background-color: rgba(0, 0, 0, 0.6);
        padding: 2rem 3rem;
        border-radius: 1rem;
        color: #f0f0f0;
        max-width: 1000px;
        margin: auto;
    }}
    .stApp {{
        background-color: #f0f2f6;
    }}
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)

style_app_container()

# -------------- Constants --------------
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "json"
STRING_METHOD = "network"

# -------------- Functions --------------
def get_string_interactions(uniprot_id, species=9606, min_score=0.4):
    params = {
        "identifiers": uniprot_id,
        "species": species,
        "caller_identity": "streamlit_app",
        "required_score": int(min_score * 1000)
    }
    url = f"{STRING_API_URL}/{STRING_OUTPUT_FORMAT}/{STRING_METHOD}"
    response = requests.post(url, data=params)
    if response.status_code == 200:
        return response.json()
    return None

def build_network(data):
    G = nx.DiGraph()
    for interaction in data:
        p1 = interaction['preferredName_A']
        p2 = interaction['preferredName_B']
        score = interaction['score']
        G.add_edge(p1, p2, weight=score)
    return G

def find_hub_genes(G, top_n=5):
    degree_dict = dict(G.degree())
    sorted_genes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)
    return [gene for gene, _ in sorted_genes[:top_n]]

def create_graph_figure(G, hub_genes):
    if not G.edges():  # Check if the graph has any edges
        return go.Figure(layout=go.Layout(title="No interactions found for the given parameters."))

    pos = nx.spring_layout(G, seed=42)
    degrees = dict(G.degree())

    edge_x, edge_y = [], []
    for src, dst in G.edges():
        x0, y0 = pos[src]
        x1, y1 = pos[dst]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='none',
        mode='lines'
    )

    node_x, node_y, node_size, node_color, node_text = [], [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        degree = degrees[node]
        node_x.append(x)
        node_y.append(y)
        node_size.append(15 + degree * 2)
        node_color.append('red' if node in hub_genes else 'royalblue')
        node_text.append(f"{node}<br>Degree: {degree}")

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[node for node in G.nodes()],
        textposition="middle center",
        textfont=dict(color='white', size=10),
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color='white')),
        hoverinfo='text',
        hovertext=node_text
    )

    layout = go.Layout(
        title="Protein Interaction Network",
        titlefont=dict(color='#333'),
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=20, r=20, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='rgba(255,255,255,1)',
        paper_bgcolor='rgba(255,255,255,1)',
    )

    return go.Figure(data=[edge_trace, node_trace], layout=layout)

# -------------- Streamlit UI --------------
st.title("🧬 Prot'n'Hub – Protein Interaction & Hub Gene Explorer")

tabs = st.tabs(["Home", "About"])
with tabs[0]:
    st.header("Explore Protein Network")

    user_input = st.text_area("Enter Protein Name or UniProt ID", height=120)

    st.subheader("Species Selection")
    species_dict = {
        "Human (Homo sapiens)": 9606,
        "Mouse (Mus musculus)": 10090,
        "Rat (Rattus norvegicus)": 10116,
        "Zebrafish (Danio rerio)": 7955,
        "Fruit fly (Drosophila melanogaster)": 7227,
        "Custom (enter manually)": None,
    }
    selected_species = st.selectbox("Choose species", list(species_dict.keys()), index=0)
    if selected_species == "Custom (enter manually)":
        species = st.number_input("Enter NCBI Taxonomy ID", value=9606)
    else:
        species = species_dict[selected_species]

    score_threshold = st.slider("Interaction Score Threshold", 0.0, 1.0, 0.4, 0.05)

    if st.button("🔍 Analyze"):
        with st.spinner("Fetching and analyzing data..."):
            uniprot_id = user_input.strip()
            if not uniprot_id:
                st.warning("Please enter a Protein Name or UniProt ID.")
                st.stop()

            string_data = get_string_interactions(uniprot_id, species, score_threshold)
            if not string_data:
                st.error("❌ No interaction data found or error fetching data. Check your input and try again.")
                st.stop()

            G = build_network(string_data)
            if G.number_of_nodes() == 0:
                st.warning("ℹ️ No network could be built with the given parameters. Try a lower score threshold or a different protein.")
                st.stop()

            hub_genes = find_hub_genes(G)
            if hub_genes:
                st.success(f"Top Hub Genes: {', '.join(hub_genes)}")
            else:
                st.info("No distinct hub genes found based on the current network.")

            fig = create_graph_figure(G, hub_genes)
            st.plotly_chart(fig, use_container_width=True)

            # Save to PNG and offer download
            buf = io.BytesIO()
            try:
                fig.write_image(buf, format="png", width=1000, height=800, engine="kaleido")
                st.download_button(
                    label="📥 Download Network as PNG",
                    data=buf.getvalue(),
                    file_name=f"{uniprot_id}_network.png",
                    mime="image/png"
                )
            except Exception as e:
                st.error(f"Error saving the figure: {e}")
                st.info("Please ensure 'kaleido' is installed. You can try: `pip install -U kaleido`")


            with st.expander("📊 Network Analysis Results", expanded=True):
                st.write(f"⭐ **Nodes**: {G.number_of_nodes()}")
                st.write(f"List of Nodes: {', '.join(list(G.nodes()))}")
                st.write(f"⭐ **Edges**: {G.number_of_edges()}")
                degree_dict = dict(G.degree())
                st.write("⭐ **Node Degrees:**")
                for node, degree in sorted(degree_dict.items(), key=lambda item: item[1], reverse=True)[:10]:
                    st.write(f"- {node}: {degree}")
                if len(degree_dict) > 10:
                    st.write("... and more.")

                if degree_dict:
                    main_hub = max(degree_dict, key=degree_dict.get)
                    st.info(f"🏆 **Main Hub Gene (Highest Degree):** {main_hub} (Degree: {degree_dict[main_hub]})")
                else:
                    st.info("No nodes to determine a main hub gene.")


with tabs[1]:
    st.header("About Prot'n'Hub")
    st.markdown("""
    ... (rest of your About section)
    """)
