import numpy as np 

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
    def __init__(self,letters, match_score, mismatch_soce, gap):
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
        distance_dict = self.distance()
        # update the line: S_i0 = i*g    
        score_m[0] = np.arange(0,-m,-1)
        # update the column: S_0j = j*g
        score_m[:,0] = np.arange(0,-n,-1)
        # update others values
        for i in range(1,n):
            for j in range(1,m):
                score_m[i,j] = max(score_m[i-1,j-1]+ distance_dict[A[j-1]][B[i-1]],
                                score_m[i-1,j  ] + self.gap,
                                score_m[i  ,j-1] + self.gap               
                )
        return score_m
    
    # plot table with direction 

    
    # Needleman_and_wunsch Algorithm wihout traceback
    def Needleman_and_wunsch_algo(self, A, B, distance_matrix, opening_gap=True):
        n = len(B)+1
        m = len(A)+1
        score_m = np.zeros((n,m))
        # This matrix contain at position index i,j the max value   
        max_drc =  np.zeros((n,m))  #max direction
    
        new_A , new_B = '',''
        
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
                indices_list = np.where(np.array(values) == max_value)[0]
                # incase gap extenstion is used, flip the array to get the direction indice of a gap
                if not opening_gap:
                    # check the previous to extend the gap
                    pass

                score_m[i,j] = max_value
                max_drc[i,j] = indices_list[0] 

        score =  score_m[i,j] 

        max_drc[0,1:] = 2*np.ones(m-1)
        max_drc[1:,0] = np.ones(n-1)

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

        
        # Retrun score_m, max_drc, sequence A, sequence B, where a gap is represented as _
        return new_A, new_B, score

