import streamlit as st
import numpy as np
from markov_chains import*



def main():
    st.title('Markov Simulation App')
    st.sidebar.title('Input Parameters')

    # Sidebar inputs
    n = st.sidebar.number_input('Enter value for n:', min_value=1, step=1, value=3)
    n_f = st.sidebar.number_input('Enter value for n_f:', min_value=1, step=1, value=3)

    st.sidebar.write("Enter preferences for pref_1 (separated by commas):")
    pref_1_input = st.sidebar.text_area("pref_1", "1,2,3\n4,5,6\n7,8,9")

    st.sidebar.write("Enter preferences for pref_2 (separated by commas):")
    pref_2_input = st.sidebar.text_area("pref_2", "1,2,3\n4,5,6\n7,8,9")

    if st.sidebar.button('Run Simulation'):
        st.write(f"Running simulation for n={n} and n_f={n_f}...")

        # Process input preferences
        try:
            pref_1 = np.array([list(map(int, row.split(','))) for row in pref_1_input.split('\n')])
            pref_2 = np.array([list(map(int, row.split(','))) for row in pref_2_input.split('\n')])
        except Exception as e:
            st.error(f"Error in preferences input format: {e}")
            return

        # Debugging: Print the preferences (optional)
        # st.write("Preferences for pref_1:")
        # st.write(pref_1)
        # st.write("Preferences for pref_2:")
        # st.write(pref_2)

        # Generate all_matchings_M
        rep = 1
        All_matchings_M = comb(n, n_f)

        # Run markov_sim function
        t, n_h_stables, M, M_m, M_f, M_e, M_e_f = markov_sim(pref_1, pref_2, All_matchings_M)

        # Display results
        st.header('Results')
        st.write("Vector t:")
        st.write(t)

        # Improved display of matrices using st.expander and st.dataframe
        matrix_dict = {
            'Uniform block probability (M)': M,
            'Ackermann Best response 1 as active side (M_m)': M_m,
            'Ackermann Best response 2 as active side (M_f)': M_f,
            'Utilitarian response, 1 as active side (M_e)': M_e,
            'Utilitarian response, 2 as active side (M_e_f)': M_e_f
        }

        for matrix_name, matrix in matrix_dict.items():
            with st.expander(matrix_name):
                # Display the selected matrix
                st.dataframe(matrix)

if __name__ == '__main__':
    main()