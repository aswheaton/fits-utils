import configparser

config = configparser.ConfigParser()
config['TELESCOPE'] = {'Size': '15'}

config['DATA SETTINGS'] = {'Standard Star A': 'bd62',
                           'Standard Star B': 'bd25',
                           'Target ID': 'm52',
                           'Bands': 'g, r, u'}

with open('config.ini', 'w') as configfile:
    config.write(configfile)
