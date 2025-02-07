import math

import pandas as pd

f_threshold = 0.01


base_df = pd.read_csv('df_gaussian-tddft-svp.csv')
acetone_df = pd.read_csv('df_gaussian-tddft-svp-acetone.csv')

df = base_df.merge(acetone_df, how='inner', on='id').rename(columns=lambda x: x.replace('_x', '_gasphase').replace('_y', '_acetone'))

df = df.drop('has_failed_gasphase', axis=1)
df = df.drop('has_failed_acetone', axis=1)

new_data = []
for _ in df.to_dict(orient='records'):

    lambda_delta = None
    f_delta = None

    vis_to_vis = None
    uv_to_vis = None
    nir_to_vis = None

    if math.isnan(_['f_max_vis_acetone']) or _['f_max_vis_acetone'] < f_threshold:
        lambda_delta = 0
        f_delta = 0
        vis_to_vis = False
        uv_to_vis = False
        nir_to_vis = False

    else:
        acetone_lambda = _['lambda_max_vis_acetone']
        acetone_f = _['f_max_vis_acetone']

        if math.isnan(_['f_max_vis_gasphase']) or _['f_max_vis_gasphase'] < f_threshold:

            if (not math.isnan(_['f_max_uv_gasphase']) and _['f_max_uv_gasphase'] > f_threshold) \
                and (not math.isnan(_['f_max_nir_gasphase']) and _['f_max_nir_gasphase'] > f_threshold):

                # determine distance to uv and nir values
                uv_distance = abs(acetone_lambda - _['lambda_max_uv_gasphase'])
                nir_distance = abs(acetone_lambda - _['lambda_max_nir_gasphase'])

                if uv_distance < nir_distance:
                    lambda_delta = acetone_lambda - _['lambda_max_uv_gasphase']
                    f_delta = acetone_f - _['f_max_uv_gasphase']
                    vis_to_vis = False
                    uv_to_vis = True
                    nir_to_vis = False
                else:
                    lambda_delta = acetone_lambda - _['lambda_max_nir_gasphase']
                    f_delta = acetone_f - _['f_max_nir_gasphase']
                    vis_to_vis = False
                    uv_to_vis = False
                    nir_to_vis = True

            elif (not math.isnan(_['f_max_uv_gasphase']) and _['f_max_uv_gasphase'] > f_threshold) \
                and (math.isnan(_['f_max_nir_gasphase']) or _['f_max_nir_gasphase'] < f_threshold):

                lambda_delta = acetone_lambda - _['lambda_max_uv_gasphase']
                f_delta = acetone_f - _['f_max_uv_gasphase']
                vis_to_vis = False
                uv_to_vis = True
                nir_to_vis = False

            elif (math.isnan(_['f_max_uv_gasphase']) or _['f_max_uv_gasphase'] < f_threshold) \
                and (not math.isnan(_['f_max_nir_gasphase']) and _['f_max_nir_gasphase'] > f_threshold):

                lambda_delta = acetone_lambda - _['lambda_max_nir_gasphase']
                f_delta = acetone_f - _['f_max_nir_gasphase']
                vis_to_vis = False
                uv_to_vis = False
                nir_to_vis = True

            else:
                print(_['id'])
                continue

        else:
            lambda_delta = acetone_lambda - _['lambda_max_vis_gasphase']
            f_delta = acetone_f - _['f_max_vis_gasphase']
            vis_to_vis = True
            uv_to_vis = False
            nir_to_vis = False


    _['lambda_delta'] = lambda_delta
    _['f_delta'] = f_delta
    _['vis_to_vis'] = vis_to_vis
    _['uv_to_vis'] = uv_to_vis
    _['nir_to_vis'] = nir_to_vis

    if lambda_delta > 0:
        _['bathochromic'] = True
    else:
        _['bathochromic'] = False

    if lambda_delta < 0:
        _['hypsochromic'] = True
    else:
        _['hypsochromic'] = False

    if f_delta > 0:
        _['hyperchromic'] = True
    else:
        _['hyperchromic'] = False

    if f_delta < 0:
        _['hypochromic'] = True
    else:
        _['hypochromic'] = False

    new_data.append(_)

df = pd.DataFrame(new_data)
df.to_csv('tmQMg*.csv', index=False)

