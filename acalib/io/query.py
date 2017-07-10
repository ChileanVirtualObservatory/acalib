from astroquery.skyview import SkyView
from astropy.table import Table


def survey_table():
    dd = SkyView.survey_dict
    survey_lst = []
    type_lst = []
    for (key, val) in dd.items():
        survey_lst.extend(val)
        type_lst.extend(len(val) * [key])
    tt = Table()
    tt['survey'] = survey_lst
    tt['type'] = type_lst
    #tt.show_in_notebook()
    return tt

def query_surveys(position,survey):
    return(SkyView.get_image_list(position=position,survey=survey))