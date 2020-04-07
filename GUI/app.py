import streamlit as st
from openSIR_GUI_test import displayApp

dataSource = st.sidebar.selectbox("Which data should be used?",
                                  ('', 'Default data', 'User provided'))

displayApp(dataSource)
