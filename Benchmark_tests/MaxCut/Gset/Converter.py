import numpy as np

# return JMatrix, HVector, nodes
def create_JMatrix_and_HVector(nome_file):
    try:
        fileToConvert = open(nome_file, 'r')
    except:
        print("Error: the file " + nome_file + " doesn't exist or cannot be open\n")
        return np.zeros((0, 0)),  np.zeros(0), 0, 0 
        
    first = True

    offset = 0

    for line in fileToConvert:	
        numbers = line.rstrip().split(" ")
        if first:
            first = False	
            if len(numbers) != 2:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
            try: 
                nodes = int(numbers[0])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
            JMatrix = np.zeros((nodes, nodes))
            HVector = np.zeros(nodes)
        else:
            if len(numbers) != 3:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
            try: 
                i = int(numbers[0])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
                
            try: 
                j = int(numbers[1])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
                
            try: 
                w = float(numbers[2])/2
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
                
            if i > nodes or j > nodes or i < 1 or j < 1:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)),  np.zeros(0), 0, 0
            
            JMatrix.itemset((i-1, j-1), w)
            offset -= w
			
    JMatrix = (JMatrix + JMatrix.T)/2
	
    return JMatrix, HVector, offset, nodes


# return QMatrix, nodes
def create_QMatrix(nome_file):
    try:
        fileToConvert = open(nome_file, 'r')
    except:
        print("Error: the file " + nome_file + " doesn't exist or cannot be open\n")
        return np.zeros((0, 0)), 0

    first = True

    offset = 0

    for line in fileToConvert:
        numbers = line.rstrip().split(" ")
        if first:
            first = False
            if len(numbers) != 2:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0
            try:
                nodes = int(numbers[0])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0
            QMatrix = np.zeros((nodes, nodes))
        else:
            if len(numbers) != 3:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0
            try:
                i = int(numbers[0])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0

            try:
                j = int(numbers[1])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0

            try:
                w = float(numbers[2])
            except:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0

            if i > nodes or j > nodes or i < 1 or j < 1:
                print("Error: file format not correct\n")
                print(line)
                print(numbers)
                return np.zeros((0, 0)), 0

            QMatrix.itemset((i - 1, j - 1), 2*w)
            QMatrix[i - 1, i - 1] -= w
            QMatrix[j - 1, j - 1] -= w

    return QMatrix, nodes
    
def move_J_matrix_to_dict(JMatrix, HVector, nodes):
    J = {}
    h = {}
    
    for i in range(nodes):
        for j in range(nodes):
            J[(i,j)] = JMatrix.item((i,j))
        h[i] = HVector.item(i)

    return J, h
