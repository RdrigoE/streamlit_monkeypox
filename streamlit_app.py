import base64
from matplotlib import pyplot as plt
import streamlit as st
import pandas as pd
from bar_plot import get_bar_chart, scatter
from get_info_structure import get_info_from_structure

apobec = st.sidebar.selectbox(
    "Select the software:",
    options=["snippy", "iVar"],
)
if apobec == "snippy":
    df = pd.read_csv("./software/snippy/monkeypox_data.csv")
else:
    df = pd.read_csv("./software/ivar/monkeypox_data.csv")

st.sidebar.header("Please filter Here:")
dataset = st.sidebar.multiselect(
    "Select the Datasets:", options=df["dataset"].unique(), default=df["dataset"].unique()
)
type = st.sidebar.multiselect(
    "Select the Type:", options=df["type"].unique(), default=df["type"].unique()
)

annotation = st.sidebar.multiselect(
    "Select the Annotation:",
    options=df["mutation_anotation"].unique(),
    default=df["mutation_anotation"].unique()
)
apobec = st.sidebar.multiselect(
    "Select the Apobec:",
    options=df["apobec"].unique(),
    default=df["apobec"].unique()
)
bins = st.sidebar.slider('Number of bins', 0, 1000, 50)
default_value = 0.0
min_freq_value = st.sidebar.slider('Min Freq', 0.0, 1.0, value=0.0)
max_freq_value = st.sidebar.slider('Max Freq', 0.0, 1.0,  value=1.0)
min_depth_value = st.sidebar.slider('Min Depth', 0, int(df["depth"].max()))
limit = st.sidebar.slider('Limit', 1.0, 0.0)

tab1, tab2, tab3, tab4, tab5 = st.tabs(
    ["plots",
     "information",
     "snippy tree",
     "ivar tree",
     "chi-test"])
with tab1:

    df_selection = df.query(
        "dataset == @dataset & depth > @min_depth_value & type == @type & mutation_anotation == @annotation & apobec == @apobec & freq >= @min_freq_value & freq <= @max_freq_value "

    )

    fig, positions = get_bar_chart(
        df_selection["pos"].value_counts(), int(bins), "Monkeypox SNPs Frequency", "monkeypox.png")
    st.pyplot(fig)
    st.info("Positions range with the most occurrences of snps")
    st.write(positions)
    st.text(f"The number of total rows is {df.shape[0]}")
    st.text(f"The number of current rows is {df_selection.shape[0]}")


# Filter the rows where the 'freq' column is greater than or equal to 0.50
    filtered_df = df_selection[df_selection['freq'] >= 0.51]

# Group the rows by the 'sample_name' column and get the number of entries for each sample
    entries_per_sample = filtered_df.groupby('sample_name').size()

# Calculate the mean of the number of entries per sample
    mean_entries = entries_per_sample.mean()

    st.text(
        f"The number of mean validated variants per sample is {round(mean_entries, 2)}")

    filtered_df = df_selection[df_selection['freq'] < 0.51]

# Group the rows by the 'sample_name' column and get the number of entries for each sample
    entries_per_sample = filtered_df.groupby('sample_name').size()

# Calculate the mean of the number of entries per sample
    mean_entries = entries_per_sample.mean()

    st.text(
        f"The number of mean minor validated variants per sample is {round(mean_entries, 2)}")
    freq_axis = df_selection["freq"]
    depth_axis = df_selection["depth"]
    fig = scatter(depth_axis, freq_axis, "Freq vs depth", "Depth", "Frequency")

    # fig = scatter(freq_axis, depth_axis,"Freq vs depth")
    st.pyplot(fig)

    filtered = df_selection.groupby(['pos', 'alt', 'ref']).size().to_dict()

    all_samples = df_selection.groupby('sample_name').size().count()
    pos = [x[0] for x in filtered.keys()]
    values = [value/all_samples for value in filtered.values()]
    fig, ax = plt.subplots()

    st.title("Percentage of samples with variants in determined positions")
    filtered_apobec = df_selection.groupby(
        ['pos', 'alt', 'ref', 'apobec']).size().to_dict()

    pos_apobec = [x[0] for x in filtered_apobec.keys()]
    values_apobec = [
        value/all_samples for value in filtered_apobec.values()]

    ax.scatter(pos, values, marker="o", edgecolors="face",
               color="#A4BF90", s=5, )
    ax.scatter(pos_apobec, values_apobec, marker="o", edgecolors="face",
               color="#3686C9", s=5, )
    ax.set_ylim(-0.1, 1.1)
    ax.axhline(1, linewidth=0.1)
    ax.axhline(0, linewidth=0.1)
    ax.axhline(limit, linewidth=0.5, color="red")

    ax.legend(["Not Apobec", "Apobec"])
    st.pyplot(fig)

    count = 0
    for v in values:
        if v >= limit:
            count += 1
    st.text(
        f"There are {count} variants that happen in at least than {limit * 100}% samples")

    # Y => Number of variants
    # X => Average depth per sample (Sum depth / number of positions)

    result = df_selection.groupby('sample_name')['depth'].sum(
    )/df_selection.groupby("sample_name").size()
    fig, ax = plt.subplots()
    ax.hist(result, bins=50, color="#3686C9")
    st.title("Histogram of mean sum of depth in each sample")
    st.pyplot(fig)

    result_val = df_selection["apobec"].notnull().sum()
    result_null = df_selection["apobec"].isnull().sum()
    fig, ax = plt.subplots()

    ax.bar(["APOBEC", "NOT APOBEC"], [result_val,
           result_null], color=["#3686C9", "#A4BF90"])
    st.title("Presence of APOBEC mutations")
    st.pyplot(fig)

    result = df_selection.groupby("mutation_anotation").size().to_dict()
    fig, ax = plt.subplots()
    ax.bar(result.keys(), result.values())
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    st.title("Annotation")
    st.pyplot(fig)

    st.dataframe(df_selection)

    result = df_selection.set_index("sample_name")['sample_depth'].to_dict()

    fig, ax = plt.subplots()
    ax.hist(result.values(), bins=100)
    ax.set_xlabel("depth")
    # for tick in ax.get_xticklabels():
    #     tick.set_rotation(90)
    st.title("Sample Mean Depth Histogram")
    st.pyplot(fig)

    # Scatter Depth vs N variants
    grouped_df = df_selection.groupby(['sample_name', 'sample_depth']
                                      ).size().reset_index(name='count')
    result = grouped_df.set_index('sample_depth')['count'].to_dict()
    fig = scatter(result.keys(), result.values(),
                  "sample_depth vs number of variants", "Depth", "Number of SNPS")
    st.title("Number of Variants vs Sample Mean Depth")
    st.pyplot(fig)


with tab2:
    info_dict = get_info_from_structure(df_selection)
    st.write(info_dict)

with tab3:
    st.title("Original Snippy Tree")
    with open("./software/snippy/validated_tree.pdf", "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)

    st.title("All Snippy Minor Variants Tree")
    with open("./software/snippy/minor_tree.tree.pdf", "rb") as f_minor:
        base64_pdf_minor = base64.b64encode(f_minor.read()).decode('utf-8')
        pdf_display_minor = F'<iframe src="data:application/pdf;base64,{base64_pdf_minor}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display_minor, unsafe_allow_html=True)

    st.title("All Snippy Minor Variants Tree with Freq above 0.10")
    with open("./software/snippy/tree_with_minor_10_snippy.tree.pdf", "rb") as f_minor:
        base64_pdf_minor = base64.b64encode(f_minor.read()).decode('utf-8')
        pdf_display_minor = F'<iframe src="data:application/pdf;base64,{base64_pdf_minor}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display_minor, unsafe_allow_html=True)


with tab4:
    st.title("Original iVar Tree")
    with open("./software/ivar/ivar.tree.pdf", "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)

    st.title("All iVar Minor Variants Tree")
    with open("./software/ivar/ivar_minor.tree.pdf", "rb") as f_minor:
        base64_pdf_minor = base64.b64encode(f_minor.read()).decode('utf-8')
        pdf_display_minor = F'<iframe src="data:application/pdf;base64,{base64_pdf_minor}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display_minor, unsafe_allow_html=True)

    st.title("All iVar Minor Variants Tree with Freq above 0.10")
    with open("./software/ivar/tree_with_minor_10_ivar.tree.pdf", "rb") as f_minor:
        base64_pdf_minor = base64.b64encode(f_minor.read()).decode('utf-8')
        pdf_display_minor = F'<iframe src="data:application/pdf;base64,{base64_pdf_minor}" width="700" height="1000" type="application/pdf"></iframe>'
        st.markdown(pdf_display_minor, unsafe_allow_html=True)
with tab5:
    table_md = '''
|   Snippy/iVar   | True | False |
|---------|------|-------|
|  **True**   | 20867|  269  |
|  **False** |  35  |   0   |
'''

# Display the table using st.markdown()
    st.markdown(table_md)

    st.markdown(
        """
Chi-squared: 0.0000
p-value: 1.0000
            """
    )
