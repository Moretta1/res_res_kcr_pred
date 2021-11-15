from subscript import *
'''
input: residue-residue contact prediction result
each column meaning:
1st number = residue 1
2nd number = residue 2
3rd and 4th number = d1 and d2 from CASP format (legacy format requirement. not relevant to predictions here)
5th number = your prediction! A probability (between 0 and 1) that the residues are in contact. With 0 = no contact and 1 = contact

output: .txt file of residue-residue pair composition(RRPC)
'''
col_names = ["label", "AA", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AK", "AL", "AM", "AN", "AP", "AQ", "AR", "AS",
             "AT", "AV", "AW", "AY", "CA", "CC", "CD", "CE", "CF", "CG", "CH", "CI", "CK", "CL", "CM", "CN", "CP",
             "CQ", "CR", "CS", "CT", "CV", "CW", "CY", "DA", "DC", "DD", "DE", "DF", "DG", "DH", "DI", "DK", "DL",
             "DM", "DN", "DP", "DQ", "DR", "DS", "DT", "DV", "DW", "DY", "EA", "EC", "ED", "EE", "EF", "EG", "EH",
             "EI", "EK", "EL", "EM", "EN", "EP", "EQ", "ER", "ES", "ET", "EV", "EW", "EY", "FA", "FC", "FD", "FE",
             "FF", "FG", "FH", "FI", "FK", "FL", "FM", "FN", "FP", "FQ", "FR", "FS", "FT", "FV", "FW", "FY", "GA",
             "GC", "GD", "GE", "GF", "GG", "GH", "GI", "GK", "GL", "GM", "GN", "GP", "GQ", "GR", "GS", "GT", "GV",
             "GW", "GY", "HA", "HC", "HD", "HE", "HF", "HG", "HH", "HI", "HK", "HL", "HM", "HN", "HP", "HQ", "HR",
             "HS", "HT", "HV", "HW", "HY", "IA", "IC", "ID", "IE", "IF", "IG", "IH", "II", "IK", "IL", "IM", "IN",
             "IP", "IQ", "IR", "IS", "IT", "IV", "IW", "IY", "KA", "KC", "KD", "KE", "KF", "KG", "KH", "KI", "KK",
             "KL", "KM", "KN", "KP", "KQ", "KR", "KS", "KT", "KV", "KW", "KY", "LA", "LC", "LD", "LE", "LF", "LG",
             "LH", "LI", "LK", "LL", "LM", "LN", "LP", "LQ", "LR", "LS", "LT", "LV", "LW", "LY", "MA", "MC", "MD",
             "ME", "MF", "MG", "MH", "MI", "MK", "ML", "MM", "MN", "MP", "MQ", "MR", "MS", "MT", "MV", "MW", "MY",
             "NA", "NC", "ND", "NE", "NF", "NG", "NH", "NI", "NK", "NL", "NM", "NN", "NP", "NQ", "NR", "NS", "NT",
             "NV", "NW", "NY", "PA", "PC", "PD", "PE", "PF", "PG", "PH", "PI", "PK", "PL", "PM", "PN", "PP", "PQ",
             "PR", "PS", "PT", "PV", "PW", "PY", "QA", "QC", "QD", "QE", "QF", "QG", "QH", "QI", "QK", "QL", "QM",
             "QN", "QP", "QQ", "QR", "QS", "QT", "QV", "QW", "QY", "RA", "RC", "RD", "RE", "RF", "RG", "RH", "RI",
             "RK", "RL", "RM", "RN", "RP", "RQ", "RR", "RS", "RT", "RV", "RW", "RY", "SA", "SC", "SD", "SE", "SF",
             "SG", "SH", "SI", "SK", "SL", "SM", "SN", "SP", "SQ", "SR", "SS", "ST", "SV", "SW", "SY", "TA", "TC",
             "TD", "TE", "TF", "TG", "TH", "TI", "TK", "TL", "TM", "TN", "TP", "TQ", "TR", "TS", "TT", "TV", "TW",
             "TY", "VA", "VC", "VD", "VE", "VF", "VG", "VH", "VI", "VK", "VL", "VM", "VN", "VP", "VQ", "VR", "VS",
             "VT", "VV", "VW", "VY", "WA", "WC", "WD", "WE", "WF", "WG", "WH", "WI", "WK", "WL", "WM", "WN", "WP",
             "WQ", "WR", "WS", "WT", "WV", "WW", "WY", "YA", "YC", "YD", "YE", "YF", "YG", "YH", "YI", "YK", "YL",
             "YM", "YN", "YP", "YQ", "YR", "YS", "YT", "YV", "YW", "YY"]

file = "MP007031_results"
win_size = 10
path = file + "/cont.txt"
content = pd.read_table(path, header=0)
content = np.array(content)[1:, :]

segments = list()
residue, marked_line = findMarkedLine(file, win_size)
b1 = list()

for i in range(0, len(marked_line)):

    # for i in range(0, 1):
    lower = marked_line[i, 3]
    upper = marked_line[i, 4]
    seg = marked_line[i, 0]
    if type(seg) is not str and isnan(seg):
        print('nan')
        # continue
    ty = marked_line[i, 1]  # positive or negative

    tmp = pd.DataFrame(np.zeros((1, 401)))

    tmp.columns = col_names

    for j in range(0, len(content)):
        prob = float(str.split(arr[4], sep=' ')[0])

        if prob < 0.8:  # threshold here: 0.8
            continue

        tmp["label"] = ty
        arr = str.split(str(content[j, 0]), sep=' ')

        # residue 1 and residue 2
        index1 = int(arr[0]) - 1
        index2 = int(arr[1]) - 1

        if index1 > upper or index2 < lower or residue[index1] is '-' or residue[index2] is '-':
            continue

        pair = str(residue[index1]) + str(residue[index2])
        print(pair)
        original_num = tmp[pair]
        tmp[pair] = original_num + 1

    test = np.asarray(tmp)
    b1.append(test)

b1 = np.asarray(b1)
b1 = b1.reshape(-1, 401)
df = pd.DataFrame(b1, columns=col_names)

# save the rrpc file（contact files）
output = "rrpc_10.txt"
df.to_csv(output, index=True)
