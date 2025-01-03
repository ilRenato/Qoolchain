def DIMACStoList(filename):
    try:
        f = open(filename, "r")
    except:
        print(f'File {filename} not found')
        exit(1)

    graph = []
    n_nodes = 0
    n_edges = 0
    for line in f:
        if line.strip() == '':
            continue
        if line[0] == 'c':
            continue
        else:
            splitted_line = line.split()
            if splitted_line[0] == 'p':
                n_nodes = int(splitted_line[2])
                n_edges = int(splitted_line[3])
            elif splitted_line[0] == 'e':
                graph.append((int(splitted_line[1]), int(splitted_line[2])))
            else:
                continue

    f.close()
    return n_nodes, n_edges, graph
