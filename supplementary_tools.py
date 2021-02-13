from datetime import datetime
import datetime as dt

def get_init_time(model):

    current_time = datetime.utcnow()
    year = current_time.year
    month = current_time.month
    day = current_time.day
    hour = current_time.hour

    if model=='HRRR':
        if hour <3:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<9:
            init_hour = '00'
        elif hour<15:
            init_hour = '06'
        elif hour<21:
            init_hour = '12'
        else:
            init_hour = '18'

    elif model=='NAM':
        if hour <4:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<10:
            init_hour = '00'
        elif hour<16:
            init_hour = '06'
        elif hour<22:
            init_hour = '12'
        else:
            init_hour = '18'

    elif model=='GFS':
        if hour <5:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<11:
            init_hour = '00'
        elif hour<17:
            init_hour = '06'
        elif hour<23:
            init_hour = '12'
        else:
            init_hour = '18'

    mdate = str(year)+str(month)+str(day)
    output = [mdate,init_hour]

    return output

def wet_bulb(temp,dewpoint):
    tdd = temp-dewpoint
    wet_bulb = temp-((1/3)*tdd)
    return wet_bulb

def fram(ice,wet_bulb,velocity):
    ilr_p = ice
    ilr_t = (-0.0071*(wet_bulb**3))-(0.039*(wet_bulb**2))-(0.3904*wet_bulb)+0.5545
    ilr_v = (0.0014*(velocity**2))+(0.0027*velocity)+0.7574

    cond_1 = np.ma.masked_where(wet_bulb>-0.35,ice)
    cond_2 = np.ma.masked_where((wet_bulb<-0.35) & (velocity>12.),ice)
    cond_3 = np.ma.masked_where((wet_bulb<-0.35) & (velocity<=12.),ice)

    cond_1 = cond_1.filled(0)
    cond_2 = cond_2.filled(0)
    cond_3 = cond_3.filled(0)

    ilr_1 = (0.7*ilr_p)+(0.29*ilr_t)+(0.01*ilr_v)
    ilr_2 = (0.73*ilr_p)+(0.01*ilr_t)+(0.26*ilr_v)
    ilr_3 = (0.79*ilr_p)+(0.2*ilr_t)+(0.01*ilr_v)

    accretion_1 = cond_1*ilr_1
    accretion_2 = cond_2*ilr_2
    accretion_3 = cond_3*ilr_3

    total_accretion=accretion_1+accretion_2+accretion_3
    return total_accretion

def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
