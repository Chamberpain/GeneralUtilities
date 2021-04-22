import os
from data_save_utilities.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead
from data_save_utilities.pickle_utilities import load,save
import datetime
import geopy
import re
import numpy as np 
import geopy.distance

class MetaClass():
    """ class to organize all of read in meta net cdf data
        ----------
        file_path: the file path of the net cdf file you wish to read

    """
    def __init__(self,profile_list):
        self.id = profile_list[0].group('id')


class ProfClass():
    """ class to organize all of read in profile net cdf data
        ----------
        file_path: the file path of the net cdf file you wish to read

    """     
    def __init__(self,profile_list):

        mask = [True]*len(profile_list)
        self.date = self.Date(profile_list,mask)
        self.pos = self.Position(profile_list,mask)
        self.speed = Speed(self.date,self.pos,speed_limit=10)
        self.depth = self.Depth(profile_list,mask)

    class Date(object):
        def __init__(self,profile_list,mask):
            self._mask = mask
            self._list = [self.parse_time(profile) for profile in profile_list]

        def parse_time(self,profile):
            """
            Function parses a netcdf time related instance and returns a time instance

            Parameters
            ----------
            variable instance of the time you wish decyfered
            reference date for juld 

            Returns
            -------
            date time instance
            """
            fmt = '%b %d %Y %H:%M:%S'
            date_string = "{} {} {} {}".format(profile.group('month'),profile.group('day'),profile.group('year'),profile.group('time'))
            return datetime.datetime.strptime(date_string,fmt)




        def return_mask(self):
            return self._mask

        def is_problem(self):
            """
            Returns boolean value to show that all the date tests have been passed.

            Current tests:
            All time differences are greater than 0
            All profiles happen before today 
            All profiles happen after the beginning of the argo program
            The maximum time difference between profiles cannot be greater than 9 months

            Returns
            -------
            Boolean value to show if this profile is a problem
            """

            time_diff_list = [self._list[idx+1]-self._list[idx] for idx in range(len(self._list)-1)]
            seconds_diff_list = [_.days*24+_.seconds/3600. for _ in time_diff_list]
            test_1 = (np.array(seconds_diff_list)>0).all()  # make sure all time is going in the right direction
            test_2 = ((np.array(self._list))<datetime.datetime.today()).all() # make sure all profiles happen before today
            test_3 = ((np.array(self._list))>datetime.datetime(2010,1,1)).all() # make sure all profiles are after beginning of argo program
            test_4 = ((np.array([_.days for _ in time_diff_list])<270)).all() # make sure maximum time difference between profiles is less than 9 months
            bool_value = ~(test_1&test_2&test_3&test_4)
            if bool_value: 
                print 'date is the problem'
            return bool_value


    class Position(BasePosition):
        """ class that is a holder of the position information

        """     
        def __init__(self,profile_list,mask):         
            self._mask = mask
            self._list = [geopy.Point(_.group('lat'),_.group('lon')) for _ in profile_list]



    class Depth():
        def __init__(self,profile_list,mask):
            self._mask = mask
            self._list = [_.group('depth') for _ in profile_list]


class WetLabsRead(BaseRead):
    data_description = 'wetlabs'
    def __init__(self,file_path):
        open_file = open(file_path,'r')
        kml_text = open_file.read() 
        open_file.close()
        split_text = re.split('</description>|<description>',kml_text)
        match_list = [re.match('NAVIS Float (?P<id>\w+) Profile (?P<profile>\w+)<br/>(?P<month>\S+) (?P<day>\S+) (?P<year>\S+) (?P<time>\S+)<br/>LAT   (?P<lat>\S+) LONG (?P<lon>\S+)<br/>DEPTH (?P<depth>\S+) m',_) 
        for _ in split_text]
        match_list = [_ for _ in match_list if _]        
        self.prof = ProfClass(match_list)
        self.meta = MetaClass(match_list)
        super(WetLabsRead,self).__init__()


def aggregate_wetlabs_list(num=-1,list_save=False,recompile=False):
    """
    function that returns dictionary of argo read classes

    Parameters
    ----------
    num: the length of the dictionary you want (default is all files)
    list_save: whether you want to save the list (default is false)

    Returns
    -------
    dictionary of argo read classes.
    """
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    def compile_matches():
        base_folder = '/Users/pchamberlain/iCloud/Data/Raw/Seabird/'
        matches = []
        for root, dirnames, filenames in os.walk(base_folder):
            for filename in filenames:
                if filename.endswith('.kml'):
                    matches.append(os.path.join(root,filename))
                print len(matches)
        save("seabird_match.pkl",matches)
        return matches

    # depth = Depth()
    if recompile:
        compile_matches()
    else:
        try:
            matches = load(__location__+'/seabird_match.pkl')
        except IOError:
            matches = compile_matches()

    num = 0
    while num<len(matches):
        print matches[num]
        yield WetLabsRead(matches[num])
        num +=1 