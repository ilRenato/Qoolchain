def evaluateResult(f, res_cpp, off):
    cpp_value = off
    for i in range(len(f)):
        for j in range(len(f)):
            if res_cpp[i] and res_cpp[j]:
                cpp_value += f[i][j]
    return cpp_value
