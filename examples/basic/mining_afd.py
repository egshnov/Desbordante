import desbordante

TABLE = 'datasets/inventory_afd.csv'
ERROR = 0.1

algo = desbordante.afd.algorithms.Tane()
algo.load_data(table=(TABLE, ',', True))
algo.execute(error=ERROR)
result = algo.get_fds()
print('AFDs:')
for fd in result:
    print(fd)
