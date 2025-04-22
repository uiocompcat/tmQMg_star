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
    vis_to_uv = None
    vis_to_nir = None

    if not math.isnan(_['f_max_vis_gasphase']) and _['f_max_vis_gasphase'] > f_threshold \
       and not math.isnan(_['f_max_vis_acetone']) and _['f_max_vis_acetone'] > f_threshold:

        print('Shift in the visible')

        lambda_delta = _['lambda_max_vis_acetone'] - _['lambda_max_vis_gasphase']
        f_delta = _['f_max_vis_acetone'] - _['f_max_vis_gasphase']
        vis_to_vis = True
        uv_to_vis = False
        nir_to_vis = False
        vis_to_uv = False
        vis_to_nir = False

    elif (math.isnan(_['f_max_vis_gasphase']) or _['f_max_vis_gasphase'] <= f_threshold) \
       and (math.isnan(_['f_max_vis_acetone']) or _['f_max_vis_acetone'] <= f_threshold):

        print('No shift')

        lambda_delta = 0
        f_delta = 0
        vis_to_vis = False
        uv_to_vis = False
        nir_to_vis = False
        vis_to_uv = False
        vis_to_nir = False

    elif (math.isnan(_['f_max_vis_gasphase']) or _['f_max_vis_gasphase'] <= f_threshold) \
       and not math.isnan(_['f_max_vis_acetone'] and _['f_max_vis_acetone'] > f_threshold):

        print("Shift to the visible")

        uv_distance = None
        nir_distance = None

        if not math.isnan(_['f_max_uv_gasphase']) and _['f_max_uv_gasphase'] > f_threshold:
            uv_distance = abs(_['lambda_max_vis_acetone'] - _['lambda_max_uv_gasphase'])
        if not math.isnan(_['f_max_nir_gasphase']) and _['f_max_nir_gasphase'] > f_threshold:
            uv_distance = abs(_['lambda_max_vis_acetone'] - _['lambda_max_nir_gasphase'])

        if uv_distance is not None and nir_distance is not None:

            if uv_distance < nir_distance:
                lambda_delta = _['lambda_max_vis_acetone'] - _['lambda_max_uv_gasphase']
                f_delta = _['f_max_vis_acetone'] - _['f_max_uv_gasphase']
                vis_to_vis = False
                uv_to_vis = True
                nir_to_vis = False
                vis_to_uv = False
                vis_to_nir = False
            else:
                lambda_delta = _['lambda_max_vis_acetone'] - _['lambda_max_nir_gasphase']
                f_delta = _['f_max_vis_acetone'] - _['f_max_nir_gasphase']
                vis_to_vis = False
                uv_to_vis = False
                nir_to_vis = True
                vis_to_uv = False
                vis_to_nir = False

        elif uv_distance is None and nir_distance is None:

            lambda_delta = 0
            f_delta = 0
            vis_to_vis = False
            uv_to_vis = False
            nir_to_vis = False
            vis_to_uv = False
            vis_to_nir = False

        elif uv_distance is not None and nir_distance is None:

            lambda_delta = _['lambda_max_vis_acetone'] - _['lambda_max_uv_gasphase']
            f_delta = _['f_max_vis_acetone'] - _['f_max_uv_gasphase']
            vis_to_vis = False
            uv_to_vis = True
            nir_to_vis = False
            vis_to_uv = False
            vis_to_nir = False

        elif uv_distance is None and nir_distance is not None:

            lambda_delta = _['lambda_max_vis_acetone'] - _['lambda_max_nir_gasphase']
            f_delta = _['f_max_vis_acetone'] - _['f_max_nir_gasphase']
            vis_to_vis = False
            uv_to_vis = False
            nir_to_vis = True
            vis_to_uv = False
            vis_to_nir = False


    elif not math.isnan(_['f_max_vis_gasphase']) and _['f_max_vis_gasphase'] > f_threshold \
       and (math.isnan(_['f_max_vis_acetone']) or _['f_max_vis_acetone'] <= f_threshold):

        print("Shift away from the visible")

        uv_distance = None
        nir_distance = None

        if not math.isnan(_['f_max_uv_acetone']) and _['f_max_uv_acetone'] > f_threshold:
            uv_distance = abs(_['lambda_max_vis_gasphase'] - _['lambda_max_uv_acetone'])
        if not math.isnan(_['f_max_nir_acetone']) and _['f_max_nir_acetone'] > f_threshold:
            uv_distance = abs(_['lambda_max_vis_gasphase'] - _['lambda_max_nir_acetone'])

        if uv_distance is not None and nir_distance is not None:

            if uv_distance < nir_distance:
                lambda_delta = _['lambda_max_uv_acetone'] - _['lambda_max_vis_gasphase']
                f_delta = _['f_max_uv_acetone'] - _['f_max_vis_gasphase']
                vis_to_vis = False
                uv_to_vis = False
                nir_to_vis = False
                vis_to_uv = True
                vis_to_nir = False
            else:
                lambda_delta = _['lambda_max_nir_acetone'] - _['lambda_max_vis_gasphase']
                f_delta = _['f_max_nir_acetone'] - _['f_max_vis_gasphase']
                vis_to_vis = False
                uv_to_vis = False
                nir_to_vis = False
                vis_to_uv = False
                vis_to_nir = True

        elif uv_distance is None and nir_distance is None:

            lambda_delta = 0
            f_delta = 0
            vis_to_vis = False
            uv_to_vis = False
            nir_to_vis = False
            vis_to_uv = False
            vis_to_nir = False

        elif uv_distance is not None and nir_distance is None:

            lambda_delta = _['lambda_max_uv_acetone'] - _['lambda_max_vis_gasphase']
            f_delta = _['f_max_uv_acetone'] - _['f_max_vis_gasphase']
            vis_to_vis = False
            uv_to_vis = False
            nir_to_vis = False
            vis_to_uv = True
            vis_to_nir = False

        elif uv_distance is None and nir_distance is not None:

            lambda_delta = _['lambda_max_nir_acetone'] - _['lambda_max_vis_gasphase']
            f_delta = _['f_max_nir_acetone'] - _['f_max_vis_gasphase']
            vis_to_vis = False
            uv_to_vis = False
            nir_to_vis = False
            vis_to_uv = False
            vis_to_nir = True

    else:
        exit()

    _['lambda_delta'] = lambda_delta
    _['f_delta'] = f_delta
    _['vis_to_vis'] = vis_to_vis
    _['uv_to_vis'] = uv_to_vis
    _['nir_to_vis'] = nir_to_vis
    _['vis_to_uv'] = vis_to_uv
    _['vis_to_nir'] = vis_to_nir


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
df.to_csv('../tmQMg*.csv', index=False)

