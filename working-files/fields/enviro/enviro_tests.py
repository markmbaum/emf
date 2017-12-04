from emf.fields import enviro as env

for i in range(1,4):
    fn1 = 'enviro-files/Envsmpl' + str(i) + '.o01'
    fn2 = 'enviro-files/Envsmpl' + str(i) + '.o02'
    env.compare_o01_o02(fn1, fn2, path='o0-comparison')
