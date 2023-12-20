import numpy as np 
import matplotlib.pyplot as plt


def matrix_to_dict(matrix):
    distance_dict = {}
    # 2.1. fill the dictionary with values from the matrix
    for i in range(len(distance_matrix)):
        row_letter = matrix[0,:]
        distance_dict[row_letter] = {}
        for j in range(len(distance_matrix[i])):
            col_letters = matrix[:,0]
            distance_dict[row_letter][col_letters] = distance_matrix[i][j]


class distance_matrix:
    def __init__(self, file_path):
        self.file_path = file_path

    def get_dict(self):
        with open(self.file_path, 'r') as file:
            lines = file.readlines()
            # Create the matrix
            matrix = [line.split() for line in lines if not line.startswith('#')]

class Needleman_and_wunsch:
    def __init__(self,letters, match_score, mismatch_soce, gap): #, opening_gap,extenssion_gap
        self.letters = letters
        self.match_score = match_score
        self.mismatch_soce = mismatch_soce
        self.gap = gap
 
    # Distance matrix 
    def distance(self):
        # 1. Create a distance matrix
        # 1.1. get number of letters
        n = len(self.letters)
        # 1.2. initialize a matrix of size (nxn) with values mismatch_soce values
        distance_matrix = self.mismatch_soce*np.ones((n,n)).copy()
        # 1.3. fill the diagonal with match_score, since it reprensent the match values 
        np.fill_diagonal(distance_matrix, self.match_score)

        # 2. Transforme it to a dictionary (in order to get letters while geting the matrix + using string as index)
        # 2.1. create an empty dictionary
        distance_dict = {}
        # 2.1. fill the dictionary with values from the matrix
        for i in range(len(distance_matrix)):
            row_letter = self.letters[i]
            distance_dict[row_letter] = {}
            for j in range(len(distance_matrix[i])):
                col_letters = self.letters[j]
                distance_dict[row_letter][col_letters] = distance_matrix[i][j]

        return distance_dict
    
    # Plot table instade a dict 
    def show(self):
        pass 

    # score_matrix : 
    def score_matrix(self,A,B):
        n = len(B)+1
        m = len(A)+1

        score_m = np.zeros((n,m))
        # In order to do the trace back we need to know from where we got the max
        max_drc =  np.zeros((n,m,2)) 

        max_drc[0,1:,:] = [1,0]
        max_drc[1:,0,:] = [-1,0]

        distance_matrix = self.distance()

        # update the line: S_i0 = i*g    
        score_m[0] = np.arange(0,-m,-1)
        # update the column: S_0j = j*g
        score_m[:,0] = np.arange(0,-n,-1)


        # update others values
        for i in range(1,n):
            for j in range(1,m):
                values = [score_m[i-1,j-1]+ distance_matrix[A[j-1]][B[i-1]],
                          score_m[i-1,j  ] + self.gap,
                          score_m[i  ,j-1] + self.gap               
                ]
                max_value = max(values)
                score_m[i,j] = max_value
                # direction : 
                dirction_index =  np.where(np.array(values) == max_value)[0][0]
                if dirction_index == 0:
                    max_drc[i,j] = [-1,1]
                elif dirction_index == 1:
                    max_drc[i,j] = [-1,0]
                else:
                    max_drc[i,j] = [1,0] 
                 
        return score_m, max_drc
    
    # plot table with direction 

    def plot_score_drection(self,A,B,matrix_dirc):
        print('hi')            

        # Generate an nxm matrix (replace this with your own data)
        n = len(B)+1
        m = len(A)+1
        # 2 for 2D vectors, change if you have more components

        # Create a grid
        x, y = np.meshgrid(np.arange(0, n), np.arange(0, m))

        # Extract components for x and y directions from the matrix
        u = matrix_dirc[:,:,0]
        v = matrix_dirc[:,:,1]
        print(x.shape)
        print(y.shape)
        print(u.shape)
        print(v.shape)
        # Create the plot
        plt.figure(figsize=(8, 6))
        plt.quiver(x, y, u, v)
        plt.title('TraceBack direction')
        plt.xlabel('A')
        plt.ylabel('B')
        plt.show()

    
    # Needleman_and_wunsch Algorithm with traceback
    def Needleman_and_wunsch_algo(self, A, B, distance_matrix):
        n = len(B)+1
        m = len(A)+1
        score_m = np.zeros((n,m))
        # This matrix contain at position index x,y from where we get the max value at position i,j 
        # (x,y) can be (i-1,j-1) or (i-1,j) or (i,j-1) 
        max_drc =  np.zeros((n,m)) 
        new_A , new_B = '',''
        i,j = n-1,m-1
        # ==============================================================================================
        socre_m, max_drc = self.score_matrix(A,B)
        score =  score_m[i,j] 
        # ==============================================================================================


        # Traceback 
        # starting for the last value in the matrix index (i,j)
        while not(i == 0 and j ==0 ): 
            
            if max_drc[i,j] == 0:
                # No gap same letter for both A and B
                # Sequence A
                new_A = A[j-1]+new_A
                # Sequence B
                new_B = B[i-1]+new_B
                # Update i and j 
                i,j = i-1, j-1
                
            elif max_drc[i,j] == 1:
                # gap
                # Sequence B
                new_A = '_'+new_A
                new_B = B[i-1]+new_B
                i = i-1

            else: 
                # gap
                # Sequence A
                new_B = '_'+ new_B
                # Sequence B
                new_A = A[j-1]+ new_A
                j = j-1

        print(score)
        print(max_drc)
        # Retrun score_m, max_drc, sequence A, sequence B, where a gap is represented as _
        return new_A, new_B, score
    


    def Needleman_and_wunsch_algo_gap(self, A, B, distance_matrix, opening_gap, extended_gap):
        print('hi')
        M = np.array([])
        X = np.array([])
        Y = np.array([])
        for i in range(10):
            for j in range(50):
                M[i,j] = max(
                            M[i-1,j-1] +  distance_matrix[A[j-1]][B[i-1]],
                            X[i-1,j-1] + distance_matrix[A[j-1]][B[i-1]],
                            Y[i-1,j-1] +  distance_matrix[A[j-1]][B[i-1]]
                        )
                
                X[i,j] = max(M[i-1,j] - opening_gap, X[i-1,j] - extended_gap)
                Y[i,j] = max(M[i,j-1] - opening_gap, Y[i,j-1] - extended_gap)
        pass



    # Recouces : 
    # https://www.youtube.com/watch?v=m2wa84YU9zA