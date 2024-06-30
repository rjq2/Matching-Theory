import streamlit as st
import numpy as np
from markov_chains import*

def main():
    st.title('Markov Simulation App')
    st.subheader('This app generates the Markov Chain of a descentralized bilateral matching process, asociated on a user-defined instance. The markov chain  of a decentralized proccess depends of the probability of each blocking pair given a arbitrary matching')
    st.sidebar.title('Input Parameters')

    # Sidebar inputs
    n = st.sidebar.number_input('Enter value for $n_1$: ', min_value=2, max_value=6, step=1, value=3)
    n_f = st.sidebar.number_input('Enter value for $n_2$:', min_value=2, max_value=6, step=1, value=3)


    st.sidebar.write("Enter preferences for $pref_1$ (separated by commas):")
    pref_1_input = st.sidebar.text_area("pref_1", "1,3,3\n2,2,1\n3,1,2")

    st.sidebar.write("Enter preferences for $pref_2$ (separated by commas):")
    pref_2_input = st.sidebar.text_area("pref_2", "3,2,2\n1,3,1\n2,1,3")

    if st.sidebar.button('Run Simulation'):
        st.write(f"Running simulation for $n_1$ = {n} and $n_2$ = {n_f} ... ")

        # Process input preferences
        try:
            pref_1 = np.array([list(map(int, row.split(','))) for row in pref_1_input.split('\n')])
            pref_2 = np.array([list(map(int, row.split(','))) for row in pref_2_input.split('\n')])
        except Exception as e:
            st.error(f"Error in preferences input format: {e}")
            return

        # Generate all_matchings_M
        All_matchings_M = comb(n, n_f)

        # Run markov_sim function
        t, n_h_stables, M, M_m, M_f, M_e, M_e_f = markov_sim(pref_1, pref_2, All_matchings_M)

        # Display results
        st.header('Results')
        st.write("Vector t:")
        st.write(t)

        # Display All_matchings_M
        st.write("All Matchings (for men perspective):")
        st.write(All_matchings_M)

        # Explanation box
        with st.expander("Explanation"):
            st.write("""
                **All Matchings (All_matchings_M)**: Any row is a possible matchings between W and M
                in the sets based on the input parameters `$n_1$` and `$n_2$`, each position on any row represents the index of man, the value is the index of a woman that is paired to, if the value is 0, it means that this spcific man is alone.

                **Matrix Descriptions**:
                - **M**: Represents the uniform block probability matrix (Biro 2012).
                - **$M_{A1}$**: Ackermann Best response when 1 is the active side.
                - **$M_{A2}$**: Ackermann Best response when 2 is the active side.
                - **$M_{U1}$**: Utilitarian response when 1 is the active side.
                - **$M_{U2}$**: Utilitarian response when 2 is the active side.
            """)

        # Improved display of matrices using st.expander and st.dataframe
        matrix_dict = {
            'Uniform block probability (M)': M,
            'Ackermann Best response 1 as active side ($M_{A1}$)': M_m,
            'Ackermann Best response 2 as active side ($M_{A2}$)': M_f,
            'Utilitarian response, 1 as active side ($M_{U1}$)': M_e,
            'Utilitarian response, 2 as active side ($M_{U2}$)': M_e_f
        }

        for matrix_name, matrix in matrix_dict.items():
            with st.expander(matrix_name):
                # Display the selected matrix
                st.dataframe(matrix)

if __name__ == '__main__':
    main()