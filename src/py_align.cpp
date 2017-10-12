
/*-------------------------------------------
 * AUTHOR: Alexandre YAHI
 * AFFILIATION: Columbia University Medical Center
                Department of Biomedical Informatics
 * FALL 2015-2017
 * email: alexandre.yahi@columbia.edu
 -------------------------------------------*/

#include "py_align.h"
 #include <utility>
#include <map>
#include <sstream>

using namespace std;

// std::pair<std::string, std::string>

void  dpm_init( int ** F, int L1, int L2)
{
        F[ 0 ][ 0 ] =  0 ;

        int i=0, j=0;

        for( j = 1; j < L2; j++ )
        {
                F[ 0 ][ j ] =  0 ; // Semi-global variant
        }
        for( i = 1; i < L1; i++ )
        {
                F[ i ][ 0 ] =  0 ; // Semi-global variant
        }
}

void  P_init( int ** P, int L1, int L2)
{
        P[ 0 ][ 0 ] =  0 ;

        int i=0, j=0;

        for( j = 1; j < L2; j++ )
        {
                P[ 0 ][ j ] =  -32767 ; // Semi-global variant
        }
        for( i = 1; i < L1; i++ )
        {
                P[ i ][ 0 ] =  0 ; // Semi-global variant
        }
}

void  Q_init( int ** Q, int L1, int L2)
{
        Q[ 0 ][ 0 ] =  0 ;

        int i=0, j=0;

        for( j = 1; j < L2; j++ )
        {
                Q[ 0 ][ j ] =  0 ; // Semi-global variant
        }
        for( i = 1; i < L1; i++ )
        {
                Q[ i ][ 0 ] =  -32767 ; // Semi-global variant
        }
}


// int max_3(int fU, int fD, int fL)
// {
//     int  max = 0 ;
//
//     if(fU >= fD && fU >= fL )
//     {
//         max = fU;
//     }
//     else if(fD > fL)
//     {
//         max = fD;
//     }
//     else
//     {
//         max = fL;
//     }
//     return max;
// }

int max_2(int a, int b)
{
    int max = 0;
    if(a >= b){
        max = a;
    }
    else{
        max = b;
    }
    return max;
}


int max_3(int fU, int fD, int fL)
{
    int max = 0;
    max = max_2(fD, fL) ;
    max = max_2(fU, max) ;
    return max;
}


std::string nw_align_aff(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2,
              int gap_op,                   /* gap opening penality */
              int gap_ext                   /* gap extension penalty */
            )
    {
        // // CONST AND VAR
        // static int number_of_times = 0;
        // number_of_times++;

        // cout << number_of_times << endl;

        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        dg_op, dg_ext, ig_op, ig_ext, no_gap ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        int ** F;
        int ** P;
        int ** Q;
        F = new int * [L1];
        P = new int * [L1];
        Q = new int * [L1];
        for( int i = 0; i < L1; i++ ){
            F[ i ] = new int[L2];
            P[ i ] = new int[L2];
            Q[ i ] = new int[L2];
        }

        dpm_init(F, L1, L2);
        P_init(P, L1, L2);
        Q_init(Q, L1, L2);
        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  5;   // Match
        const int  b = -4;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = 1; j < L2; j++ )
                {
                    // Take maximum because penalties are negative
                    //Compute first P the deletion gap-opening/gap-extension matrix
                    dg_op = F[ i-1 ][ j ] - gap_op;
                    dg_ext = P[i-1][j] - gap_ext;

                    P[i][j] = max_2(dg_op, dg_ext); // DELETION

                    ig_op = F[ i ][ j-1 ] - gap_op;
                    ig_ext = Q[i][j-1] - gap_ext;

                    Q[i][j] = max_2(ig_op, ig_ext); // INSERTION

                    no_gap = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH

                    F[ i ][ j ] = max_3(no_gap, P[i][j], Q[i][j]) ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        // Look for max in the last column
        temp_i = F[0][L2-1];
        for(int k = 1; k < L1; k++)
        {
            if( F[k][L2-1] > temp_i )
            {
                i_max = k;
                temp_i = F[k][L2-1];
            }
        }

        // Look for max in the last row
        temp_j = F[L1-1][0];
        for(int l = 1; l < L2; l++)
        {
            if( F[L1-1][l] > temp_j )
            {
                j_max = l;
                temp_j = F[L1-1][l];
            }
        }

        if (temp_j > temp_i) // Start with the max of the last row
        {
          i_max = L1-1;
          while(j>j_max)
          {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
          }
        }
        else                // Start with the max of the last column
        {
          j_max = L2-1;
          while(i>i_max)
          {
              seq_1_al = seq_1[i-1] + seq_1_al;
              seq_2_al = '-' + seq_2_al;
              i--;
          }
        }

        score = F[i_max][j_max];

        int mat = 0; //0 is F, 1 is P, 2 is Q - help to jump between matrices

        // TRANSVERSAL TRACEBACK
        while( (i > 0) && (j > 0) )
        {
            if(mat == 0)
            {
                if(F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    i--;
                    j--;
                }
                else if(F[i][j] == P[i][j])
                {
                    mat = 1; // Jump to matrix P for deletion
                }
                else
                {
                    mat = 2; // Jump to matrix Q for insertion
                }
            }
            else if(mat == 1) // DELETION - matrix P
            {
                if(P[i][j]==P[i-1][j] - gap_ext) //Extension, stay in P
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 1;
                }
                else // Opening - move back to F
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 0;
                }
            }
            else // INSERTION - matrix Q
            {
                if(Q[i][j] == Q[i][j-1] - gap_ext) //Extension - stay in Q
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 2;
                }
                else //Opening - move back to F
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 0;
                }
            }
        }


        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // cout << seq_1_al << endl;
        // cout << seq_2_al << endl;

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory
        for (int i = 0; i < L1; ++i)
        {
            delete [] F[i];
            delete [] P[i];
            delete [] Q[i];
        }

        delete [] F;
        delete [] P;
        delete [] Q;

        return alignment;

        // return  std::make_pair(seq_1_al, seq_2_al) ;
}

std::string nw_align_aff_param(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2,
              int gap_op,                   /* gap opening penality */
              int gap_ext,                   /* gap extension penalty */
              int match,                    /* Match score */
              int mismatch                  /* Mismatch penalty */
            )
    {
        // // CONST AND VAR
        // static int number_of_times = 0;
        // number_of_times++;

        // cout << number_of_times << endl;

        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        dg_op, dg_ext, ig_op, ig_ext, no_gap ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        int ** F;
        int ** P;
        int ** Q;
        F = new int * [L1];
        P = new int * [L1];
        Q = new int * [L1];
        for( int i = 0; i < L1; i++ ){
            F[ i ] = new int[L2];
            P[ i ] = new int[L2];
            Q[ i ] = new int[L2];
        }

        dpm_init(F, L1, L2);
        P_init(P, L1, L2);
        Q_init(Q, L1, L2);
        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  match;   // Match
        const int  b = mismatch;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = 1; j < L2; j++ )
                {
                    // Take maximum because penalties are negative
                    //Compute first P the deletion gap-opening/gap-extension matrix
                    dg_op = F[ i-1 ][ j ] - gap_op;
                    dg_ext = P[i-1][j] - gap_ext;

                    P[i][j] = max_2(dg_op, dg_ext); // DELETION

                    ig_op = F[ i ][ j-1 ] - gap_op;
                    ig_ext = Q[i][j-1] - gap_ext;

                    Q[i][j] = max_2(ig_op, ig_ext); // INSERTION

                    no_gap = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH

                    F[ i ][ j ] = max_3(no_gap, P[i][j], Q[i][j]) ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        // Look for max in the last column
        temp_i = F[0][L2-1];
        for(int k = 1; k < L1; k++)
        {
            if( F[k][L2-1] > temp_i )
            {
                i_max = k;
                temp_i = F[k][L2-1];
            }
        }

        // Look for max in the last row
        temp_j = F[L1-1][0];
        for(int l = 1; l < L2; l++)
        {
            if( F[L1-1][l] > temp_j )
            {
                j_max = l;
                temp_j = F[L1-1][l];
            }
        }

        if (temp_j > temp_i) // Start with the max of the last row
        {
          i_max = L1-1;
          while(j>j_max)
          {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
          }
        }
        else                // Start with the max of the last column
        {
          j_max = L2-1;
          while(i>i_max)
          {
              seq_1_al = seq_1[i-1] + seq_1_al;
              seq_2_al = '-' + seq_2_al;
              i--;
          }
        }

        score = F[i_max][j_max];

        int mat = 0; //0 is F, 1 is P, 2 is Q - help to jump between matrices

        // TRANSVERSAL TRACEBACK
        while( (i > 0) && (j > 0) )
        {
            if(mat == 0)
            {
                if(F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    i--;
                    j--;
                }
                else if(F[i][j] == P[i][j])
                {
                    mat = 1; // Jump to matrix P for deletion
                }
                else
                {
                    mat = 2; // Jump to matrix Q for insertion
                }
            }
            else if(mat == 1) // DELETION - matrix P
            {
                if(P[i][j]==P[i-1][j] - gap_ext) //Extension, stay in P
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 1;
                }
                else // Opening - move back to F
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 0;
                }
            }
            else // INSERTION - matrix Q
            {
                if(Q[i][j] == Q[i][j-1] - gap_ext) //Extension - stay in Q
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 2;
                }
                else //Opening - move back to F
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 0;
                }
            }
        }


        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // cout << seq_1_al << endl;
        // cout << seq_2_al << endl;

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory
        for (int i = 0; i < L1; ++i)
        {
            delete [] F[i];
            delete [] P[i];
            delete [] Q[i];
        }

        delete [] F;
        delete [] P;
        delete [] Q;

        return alignment;

        // return  std::make_pair(seq_1_al, seq_2_al) ;
}


std::string nw_align_aff_mem(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2,
              int gap_op,                   /* gap opening penality */
              int gap_ext,                   /* gap extension penalty */
              int sim,
              int terminate
            )
    {
        // // CONST AND VAR
        // static int number_of_times = 0;
        // number_of_times++;

        // cout << number_of_times << endl;

        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        dg_op, dg_ext, ig_op, ig_ext, no_gap ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        static int ** M;
        static int ** P;
        static int ** Q;

        if(sim==0){             // Only initialize if zero similarity with previous read
            M = new int * [L1];
            P = new int * [L1];
            Q = new int * [L1];
            for( int i = 0; i < L1; i++ ){
                M[ i ] = new int[L2];
                P[ i ] = new int[L2];
                Q[ i ] = new int[L2];
            }

            dpm_init(M, L1, L2);
            P_init(P, L1, L2);
            Q_init(Q, L1, L2);
        }
        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  5;   // Match
        const int  b = -4;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = sim + 1; j < L2; j++ )
                {
                    // Take maximum because penalties are negative
                    //Compute first P the deletion gap-opening/gap-extension matrix
                    dg_op = M[ i-1 ][ j ] - gap_op;
                    dg_ext = P[i-1][j] - gap_ext;

                    P[i][j] = max_2(dg_op, dg_ext); // DELETION

                    ig_op = M[ i ][ j-1 ] - gap_op;
                    ig_ext = Q[i][j-1] - gap_ext;

                    Q[i][j] = max_2(ig_op, ig_ext); // INSERTION

                    no_gap = M[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH

                    M[ i ][ j ] = max_3(no_gap, P[i][j], Q[i][j]) ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        // Look for max in the last column
        temp_i = M[0][L2-1];
        for(int k = 1; k < L1; k++)
        {
            if( M[k][L2-1] > temp_i )
            {
                i_max = k;
                temp_i = M[k][L2-1];
            }
        }

        // Look for max in the last row
        temp_j = M[L1-1][0];
        for(int l = 1; l < L2; l++)
        {
            if( M[L1-1][l] > temp_j )
            {
                j_max = l;
                temp_j = M[L1-1][l];
            }
        }

        if (temp_j > temp_i) // Start with the max of the last row
        {
          i_max = L1-1;
          while(j>j_max)
          {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
          }
        }
        else                // Start with the max of the last column
        {
          j_max = L2-1;
          while(i>i_max)
          {
              seq_1_al = seq_1[i-1] + seq_1_al;
              seq_2_al = '-' + seq_2_al;
              i--;
          }
        }

        score = M[i_max][j_max];

        int mat = 0; //0 is M, 1 is P, 2 is Q - help to jump between matrices

        // TRANSVERSAL TRACEBACK
        while( (i > 0) && (j > 0) )
        {
            if(mat == 0)
            {
                if(M[i][j] == M[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    i--;
                    j--;
                }
                else if(M[i][j] == P[i][j])
                {
                    mat = 1; // Jump to matrix P for deletion
                }
                else
                {
                    mat = 2; // Jump to matrix Q for insertion
                }
            }
            else if(mat == 1) // DELETION - matrix P
            {
                if(P[i][j]==P[i-1][j] - gap_ext) //Extension, stay in P
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 1;
                }
                else // Opening - move back to F
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 0;
                }
            }
            else // INSERTION - matrix Q
            {
                if(Q[i][j] == Q[i][j-1] - gap_ext) //Extension - stay in Q
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 2;
                }
                else //Opening - move back to F
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 0;
                }
            }
        }


        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // cout << seq_1_al << endl;
        // cout << seq_2_al << endl;

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory if terminate
        if (terminate == 1)
        {
            for (int i = 0; i < L1; ++i)
                {
                    delete [] M[i];
                    delete [] P[i];
                    delete [] Q[i];
                }

            delete [] M;
            delete [] P;
            delete [] Q;
        }

        return alignment;
}

std::string nw_free_trail(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2
            )
    {
        // CONST AND VAR

        const int  d = 8 ;                 /* gap penalty */

        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        fU, fD, fL ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        int ** F;
        F = new int * [L1];
        for( int i = 0; i < L1; i++ )  F[ i ] = new int[L2];

        F[ 0 ][ 0 ] =  0 ;

        for( j = 1; j < L2; j++ )
        {
                F[ 0 ][ j ] =  -j*d ; // No free head gap
        }
        for( i = 1; i < L1; i++ )
        {
                F[ i ][ 0 ] =  -i*d ; // No free head hap
        }

        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  5;   // Match
        const int  b = -4;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = 1; j < L2; j++ )
                {

                    fU = F[ i-1 ][ j ] - d ; // DELETION
                    fD = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH
                    fL = F[ i ][ j-1 ] - d ; // INSERTION

                    int  max = 0 ;

                    if( fU >= fD && fU >= fL )
                    {
                            max = fU ;
                    }
                    else if( fD > fL )
                    {
                            max = fD ;
                    }
                    else
                    {
                            max = fL ;
                    }

                    F[ i ][ j ] = max ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        if(L1>L2)
        {
            temp_i = F[0][L2-1];
            for(int k = 1; k < L1; k++)
            {
                if( F[k][L2-1] > temp_i )
                {
                    i_max = k;
                    temp_i = F[k][L2-1];
                }
            }
            j_max = L2-1;

            while(i>i_max)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else if(L2>L1)
        {
            temp_j = F[L1-1][0];
            for(int l = 1; l < L2; l++)
            {
                if( F[L1-1][l] > temp_j )
                {
                    j_max = l;
                    temp_j = F[L1-1][l];
                }
            }
            i_max = L1-1;
            while(j>j_max)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }
        else{
            i_max = L1-1;
            j_max = L2-1;
        }

        score = F[i_max][j_max];

        // TRANSVERSAL TRACEBACK
        while( i > 0 || j > 0 )
        {
            if(i > 0 && j > 0 && F[i][j] == (F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]))
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                i--;
                j--;
            }
            else if(i > 0 && F[i][j] == F[ i-1 ][ j ] - d )
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
            else
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }


        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory
        for (int i = 0; i < L1; ++i)
            delete [] F[i];
        delete [] F;

        return alignment;

        // return  std::make_pair(seq_1_al, seq_2_al) ;
}

std::string nw_align_mem(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2,
              int gap_penalty,            /* gap penalty */
              int sim,                       // Similarity
             int terminate)
    {
        // CONST AND VAR
        // static int number_of_times = 0;
        // number_of_times++;

        // cout << number_of_times << endl;


        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        fU, fD, fL ;
        int        i = 0, j = 0;

        static int ** F;
        // Dynamic programming matrix - INITALIZATION
        if(sim == 0){
            F = new int * [L1];
            for( int i = 0; i < L1; i++ )  F[ i ] = new int[L2];

            F[ 0 ][ 0 ] =  0 ;

            for( j = 1; j < L2; j++ )
            {
                    F[ 0 ][ j ] =  0 ; // Semi-global variant
            }
            for( i = 1; i < L1; i++ )
            {
                    F[ i ][ 0 ] =  0 ; // Semi-global variant
            }
        }

        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  5;   // Match
        const int  b = -4;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = sim + 1; j < L2; j++ )
                {

                    fU = F[ i-1 ][ j ] - gap_penalty ; // DELETION
                    fD = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH
                    fL = F[ i ][ j-1 ] - gap_penalty ; // INSERTION

                    int  max = 0 ;

                    if( fU >= fD && fU >= fL )
                    {
                            max = fU ;
                    }
                    else if( fD > fL )
                    {
                            max = fD ;
                    }
                    else
                    {
                            max = fL ;
                    }

                    F[ i ][ j ] = max ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        // Look for max in the last column
        temp_i = F[0][L2-1];
        for(int k = 1; k < L1; k++)
        {
            if( F[k][L2-1] > temp_i )
            {
                i_max = k;
                temp_i = F[k][L2-1];
            }
        }

        // Look for max in the last row
        temp_j = F[L1-1][0];
        for(int l = 1; l < L2; l++)
        {
            if( F[L1-1][l] > temp_j )
            {
                j_max = l;
                temp_j = F[L1-1][l];
            }
        }

        if (temp_j > temp_i) // Start with the max of the last row
        {
          i_max = L1-1;
          while(j>j_max)
          {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
          }
        }
        else                // Start with the max of the last column
        {
          j_max = L2-1;
          while(i>i_max)
          {
              seq_1_al = seq_1[i-1] + seq_1_al;
              seq_2_al = '-' + seq_2_al;
              i--;
          }
        }

        score = F[i_max][j_max];

        // TRANSVERSAL TRACEBACK
        while( i > 0 && j > 0 )
        {
            if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                i--;
                j--;
            }
            else if(i > 0 && F[i][j] == F[ i-1 ][ j ] - gap_penalty )
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
            else
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // cout << seq_1_al << endl;
        // cout << seq_2_al << endl;

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory if terminate
        if (terminate == 1)
        {
        for (int i = 0; i < L1; ++i)
            delete [] F[i];
        delete [] F;
        }

        return alignment;

        // return  std::make_pair(seq_1_al, seq_2_al) ;
}

std::string nw_align(                  // Needleman-Wunsch algorithm
              string     seq_1,
              string     seq_2,
              int gap_penalty
            )
    {
        // // CONST AND VAR
        // static int number_of_times = 0;
        // number_of_times++;

        // cout << number_of_times << endl;

        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        int        fU, fD, fL ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        int ** F;
        F = new int * [L1];
        for( int i = 0; i < L1; i++ )  F[ i ] = new int[L2];

        F[ 0 ][ 0 ] =  0 ;

        for( j = 1; j < L2; j++ )
        {
                F[ 0 ][ j ] =  0 ; // Semi-global variant
        }
        for( i = 1; i < L1; i++ )
        {
                F[ i ][ 0 ] =  0 ; // Semi-global variant
        }

        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const int  a =  5;   // Match
        const int  b = -4;   // Mismatch

        const int  s[ 5 ][ 5 ] = { { a, b, b, b, -2 },    /* substitution matrix */
                                   { b, a, b, b, -2 },
                                   { b, b, a, b, -2 },
                                   { b, b, b, a, -2 },
                                   {-2,-2,-2,-2, -1 }} ;


        for( i = 1; i < L1; i++ )
        {
                for( j = 1; j < L2; j++ )
                {

                    fU = F[ i-1 ][ j ] - gap_penalty ; // DELETION
                    fD = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH
                    fL = F[ i ][ j-1 ] - gap_penalty ; // INSERTION

                    int  max = 0 ;

                    if( fU >= fD && fU >= fL )
                    {
                            max = fU ;
                    }
                    else if( fD > fL )
                    {
                            max = fD ;
                    }
                    else
                    {
                            max = fL ;
                    }

                    F[ i ][ j ] = max ;

                }
        }
        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        int temp_i = 0;
        int temp_j = 0;

        int score = 0;

        i = L1-1;
        j = L2-1;

        // Look for max in the last column
        temp_i = F[0][L2-1];
        for(int k = 1; k < L1; k++)
        {
            if( F[k][L2-1] > temp_i )
            {
                i_max = k;
                temp_i = F[k][L2-1];
            }
        }

        // Look for max in the last row
        temp_j = F[L1-1][0];
        for(int l = 1; l < L2; l++)
        {
            if( F[L1-1][l] > temp_j )
            {
                j_max = l;
                temp_j = F[L1-1][l];
            }
        }

        if (temp_j > temp_i) // Start with the max of the last row
        {
          i_max = L1-1;
          while(j>j_max)
          {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
          }
        }
        else                // Start with the max of the last column
        {
          j_max = L2-1;
          while(i>i_max)
          {
              seq_1_al = seq_1[i-1] + seq_1_al;
              seq_2_al = '-' + seq_2_al;
              i--;
          }
        }

        score = F[i_max][j_max];

        // TRANSVERSAL TRACEBACK
        while( i > 0 && j > 0 )
        {
            if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                i--;
                j--;
            }
            else if(i > 0 && F[i][j] == F[ i-1 ][ j ] - gap_penalty )
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
            else
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // cout << seq_1_al << endl;
        // cout << seq_2_al << endl;

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory
        for (int i = 0; i < L1; ++i)
            delete [] F[i];
        delete [] F;

        return alignment;

        // return  std::make_pair(seq_1_al, seq_2_al) ;
}
