from astropy.table import Table

t = Table.read('plfit.fits')

for key in ['Ncloud', 'Ncloudfit']:
    string = key
    for row in t:
        string += '{:4.0f}'.format(row[key])+' & '
    print(string+'\\\\')

for key in ['M1', 'Mean5']:
    string = key+'   '
    for row in t:
        string += '{:4.1f}'.format(row[key]/1e6) + ' & '
    print(string + '\\\\')

for base in ['index']:
    for key in ['_purepl', '_sch', '_trunc']:
        string = base + key + '   '
        for row in t:
            string += '${:4.1f}'.format(row[base + key] - 1) +\
                '^{+' + '{:4.1f}'.format(row[base + '85' + key] -
                                         row[base + key]) + '}'+ \
                '_{-' + '{:4.1f}'.format(row[base + key] -
                                         row[base + '15' + key]) + '}' + \
                                         '$ & '
        print(string + '\\\\')

for base in ['Mtrun']:
    for key in ['_purepl', '_sch', '_trunc']:
        string = base + key + '   '
        for row in t:
            string += '${:4.1f}'.format(row[base + key] / 1e6) +\
                '^{+' + '{:4.1f}'.format((row[base + '85' + key] -
                                          row[base + key]) / 1e6) + '}'+ \
                '_{-' + '{:4.1f}'.format((row[base + key] -
                                          row[base + '15' + key]) / 1e6) + '}' + \
                                         '$ & '
        print(string + '\\\\')

for base in ['pkprob']:
    for key in ['_purepl', '_sch', '_trunc']:
        string = base + key + '   '
        for row in t:
            string += '${:4.2f}'.format(row[base + key]) + '$ & '
        print(string + '\\\\')
