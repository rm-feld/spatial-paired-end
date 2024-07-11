import numpy as np
from data_synth import * 

# n = 6

# validator = Seq(n = n)
# X = validator.X + np.random.normal(0, 1/100, (n, n))

# print("real data matrix:")
# print(validator.X)

# print("\nreal data with noise:")
# print(X)

# print('-' * 30)

# forbidden_edges = []

# # find reference population
# notdone = True
# k = 0
# eps = 0.01
# while notdone and k <= 100:
#     for i in range(n):
#         p = X[i][0]
#         if all([any((X[:, j] < p + eps) & (X[:, j] > p - eps)) for j in range(1, n)]):
#             notdone = False
#             _i = i
#             break
#     k += 1
#     eps += 0.001



# _p = [p]
# _pi = [_i]
# for i in range(1, n):
#     for j in range(n):
#         if p - eps < X[j][i] < p + eps:
#             _p.append(X[j][i])
#             _pi.append(j)
#             break

# print(_p, _pi)


# # iterate through to find m
# try:
#     dic = {}
#     m = np.ones(n)
#     for i in range(n):
#         pass
# except:
#     pass


def m_cost(ref, col, m, maxeps = 0.5, debug = False, f = None):
    ''' 
    calculates an error sum for assumption of m for the column.

    returns: an error, 
    '''
    if debug and f is not None:
        print(f'data col: ({col.flatten()})', file = f)
    elif debug:
        print(f'data col: ({col.flatten()})')
    eps = 0.01
    working = True
    while eps < maxeps and working:
    
        # find all possible matches 
        counts = np.asarray([(np.abs(col[i] - ref) < eps).sum() + (np.abs(col[i] - m * ref) < eps).sum()
                         for i in range(len(col))])
        if counts.sum() < len(col):
            eps *= 1.3
        else:
            working = False

    _ref = np.copy(ref)
    mref = m * _ref
    if debug and f is not None:
        print("refs:", file = f)
        print(f'\tref: {_ref}', file = f)
        print(f'\tmref: {mref}', file = f)
    elif debug:
        print("refs:")
        print(f'\tref: {_ref}')
        print(f'\tmref: {mref}')

    # reconstruct 
    err = 0
    match = np.zeros(len(col))
    num =  np.zeros(len(col))
    # iterate through min matches: if we didn't get enough, choose a random shuffle
    if working == True:
        counts = np.random.shuffle(np.arange(len(col)))
    for i in np.argsort(counts).tolist():
        j = np.argmin(np.concatenate([np.abs(_ref - col[i]), np.abs(mref - col[i])]), axis = None)
        if debug and f is not None:
            print(f'j = {j}, grabbing value {_ref[j] if j < len(col) else mref[j % len(col)]} for {col[i]}', file = f)
        elif debug:
            print(f'j = {j}, grabbing value {_ref[j] if j < len(col) else mref[j % len(col)]} for {col[i]}')
        if j >= len(col):
            err += np.abs(col[i] - mref[j % len(col)])
            match[i] = j % len(col)
            num[i] = m
            mref[j % len(col)] = -100
            _ref[j % len(col)] = -100
        else:
            err += np.abs(col[i] - _ref[j])
            match[i] = j
            num[i] = 1
            _ref[j] = -100
            mref[j] = -100
        if debug and f is not None:     
            print('updated refs:', file = f)
            print(f'\tref: {_ref}, \n\tmref:{mref}', file = f)
            print(f'error: {err}, i: {i}, j: {j}', file = f)
            print(f'matches: {match}', file = f)
            print(f'numbers: {num}', file = f)
        elif debug:
            print('updated refs:')
            print(f'\tref: {_ref}, \n\tmref:{mref}')
            print(f'error: {err}, i: {i}, j: {j}')
            print(f'matches: {match}')
            print(f'numbers: {num}')
        
    return (err, match, num)

def m_cost_validator(nruns, id = None):
    # id for textfile
    id = np.random.randint(100, 1000)

    errs = []
    accs = []

    # run through nruns simulations
    for i in range(nruns):
        n = random.randint(3, 7)
        seq = Seq(n = n)
        X = seq.X + np.random.normal(0, 1/100, (n, n))

        ref = X[:, 0].flatten()

        for i in range(1, n):
            err, match, num = m_cost(ref, X[:, i].flatten(), seq.m[i])
            errs.append(err)
            if not np.array_equal(num, seq.I[:, i].flatten()):
                # add error to accs
                accs.append(0)
                # rerun and save to file
                with open(f'./test/mcost_debug_log_{id}.txt', 'a') as f:
                    print('-' * 20 + '\n', file = f)
                    err, match, num = m_cost(ref, X[:, i].flatten(), seq.m[i], f = f, debug = True)
                    print(f'FINAL ERROR FOR ROW: {err}', file = f)
                    print(f'\tMATCHING: {match}', file = f)
                    print(f'\tNUM: {num}, ACTUAL: {seq.I[:, i].flatten()}', file = f)
                    print('-' * 20 + '\n', file = f)
            else:
                accs.append(1)
        
    _errs = np.asarray(errs).copy()
    _accs = np.asarray(accs).copy()
    with open(f'./test/mcost_debug_log_{id}.txt', 'a') as f:
        print(f'\n\n\nSUCCESS RATE: {np.sum(_accs)/len(_accs)}', file = f)


m_cost_validator(5)
    







# x = np.asarray([0.11, 0.4, 0.2, 0.7])
# y = np.asarray([0.2, 0.84, 0.2, 0.72])

# m_cost(x, y, 2)









