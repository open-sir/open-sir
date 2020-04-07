import io
from pathlib import Path
import pandas as pd
import streamlit as st
from test_from_notebook import *

@st.cache
def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()

def load_user_data_file(datafile, file_format):
    if file_format == 'CSV':
        df = pd.read_csv(datafile)
    elif file_format == 'JSON':
        df = pd.read_json(datafile)
    elif file_format == 'XLS':
        df = pd.read_excel(datafile)
    datalist = df.iloc[:, 1].to_list()
    return datalist

# dataSource = st.sidebar.selectbox("Which data should be used?",
#                                   ('', 'Default data', 'User provided'))

def displayApp(dataSource):
    if len(dataSource) == 0:
        welcome_message = read_markdown_file("frontMatter.md")
        st.markdown(welcome_message, unsafe_allow_html=True)

    elif dataSource == 'Default data':
        #UK data from notebook
        data_list = [3269, 3983, 5018, 5683, 6650,
                     8077, 9529, 11658, 14543, 17089]

        population = st.sidebar.number_input("Population",
                                             value=67886011.0,
                                             step=1.0)

        if population > 0:

            alpha = st.sidebar.slider("Reproduction rate (alpha)",
                                      min_value=0.1,
                                      max_value=5.0,
                                      value=0.95)
            beta = st.sidebar.slider("Removal rate (beta)",
                                     min_value=0.1,
                                     max_value=5.0,
                                     value=0.38)
            infected = st.sidebar.slider("Infected",
                                         min_value=0.0,
                                         max_value=float(max(data_list)),
                                         value=float(data_list[0]),
                                         step=None,
                                         format='%i',
                                         key=None)
            recovered = st.sidebar.slider("Recovered/Removed",
                                          min_value=0.0,
                                          max_value=float(max(data_list)),
                                          value=0.0,
                                          step=None,
                                          format='%i',
                                          key=None)

            susceptible = population - infected - recovered

            if susceptible < 0:
                st.write('The number of infected plus recovered/removed' +
                         'cannot exceed the population')
            else:
                example(population, infected, recovered,
                        data_list, alpha=alpha, beta=beta)
        else:
            st.write('Population must be an integer greater than zero')


    ###########################################################################
    ###########################################################################
    elif dataSource == 'User provided':
        ext_dict = {"CSV": ("csv"),
                    "JSON": ("json"),
                    "XLS": ("xls", "xlsx")}

        user_file_ext = st.sidebar.radio("What is the format of the data file?",
                                         ('CSV', 'XLS', 'JSON'))
        user_data_file = st.sidebar.file_uploader("Upload a file",
                                                  type=ext_dict[user_file_ext])
        try:
            data_list = load_user_data_file(user_data_file, user_file_ext)

            population = st.sidebar.number_input("Population",
                                                 value=67886011.0,
                                                 step=1.0)

            if population > 0:

                alpha = st.sidebar.slider("Reproduction rate (alpha)",
                                          min_value=0.1,
                                          max_value=5.0,
                                          value=0.95)
                beta = st.sidebar.slider("Removal rate (beta)",
                                         min_value=0.1,
                                         max_value=5.0,
                                         value=0.38)
                infected = st.sidebar.slider("Infected",
                                             min_value=0.0,
                                             max_value=float(max(data_list)),
                                             value=float(data_list[0]),
                                             step=None,
                                             format='%i',
                                             key=None)
                recovered = st.sidebar.slider("Recovered/Removed",
                                              min_value=0.0,
                                              max_value=float(max(data_list)),
                                              value=0.0,
                                              step=None,
                                              format='%i',
                                              key=None)

                susceptible = population - infected - recovered

                if susceptible < 0:
                    st.write('The number of infected plus recovered/removed' +
                             'cannot exceed the population')
                else:
                    example(population, infected, recovered,
                            data_list, alpha=alpha, beta=beta)
            else:
                st.write('Population must be an integer greater than zero')

        except:
            st.write('Please upload a data file')
    return
