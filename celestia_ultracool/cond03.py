"""Interpolate physical properties from the COND03 evolutionary model grids."""
# https://perso.ens-lyon.fr/isabelle.baraffe/COND03_models

import numpy as np
from scipy.interpolate import LinearNDInterpolator

t_teff_points = np.log10(
    [[0.001, 628], [0.001, 942], [0.001, 1285], [0.001, 1553], [0.001, 1747],
     [0.001, 1901], [0.001, 2004], [0.001, 2098], [0.001, 2159], [0.001, 2207],
     [0.001, 2251], [0.001, 2321], [0.001, 2400], [0.001, 2484], [0.001, 2598],
     [0.001, 2746], [0.001, 2768], [0.001, 2824], [0.001, 2853], [0.001, 2858],
     [0.001, 2833], [0.001, 2869], [0.001, 2867], [0.001, 2856], [0.005, 455],
     [0.005, 644], [0.005, 901], [0.005, 1098], [0.005, 1263], [0.005, 1413],
     [0.005, 1543], [0.005, 1668], [0.005, 1786], [0.005, 1885], [0.005, 1965],
     [0.005, 2102], [0.005, 2264], [0.005, 2452], [0.005, 2616], [0.005, 2754],
     [0.005, 2806], [0.005, 2875], [0.005, 2916], [0.005, 2921], [0.005, 2912],
     [0.005, 2946], [0.005, 2974], [0.005, 2983], [0.01, 394], [0.01, 540],
     [0.01, 746], [0.01, 920], [0.01, 1064], [0.01, 1193], [0.01, 1316],
     [0.01, 1426], [0.01, 1536], [0.01, 1635], [0.01, 1731], [0.01, 1915],
     [0.01, 2188], [0.01, 2442], [0.01, 2622], [0.01, 2734], [0.01, 2821],
     [0.01, 2871], [0.01, 2919], [0.01, 2927], [0.01, 2942], [0.01, 2958],
     [0.01, 2993], [0.01, 3024], [0.05, 285], [0.05, 375], [0.05, 491],
     [0.05, 585], [0.05, 676], [0.05, 756], [0.05, 840], [0.05, 928],
     [0.05, 1010], [0.05, 1085], [0.05, 1171], [0.05, 1520], [0.05, 1903],
     [0.05, 1909], [0.05, 2246], [0.05, 2480], [0.05, 2651], [0.05, 2759],
     [0.05, 2837], [0.05, 2852], [0.05, 2872], [0.05, 2901], [0.05, 2948],
     [0.05, 2989], [0.1, 240], [0.1, 309], [0.1, 425], [0.1, 493], [0.1, 563],
     [0.1, 630], [0.1, 688], [0.1, 760], [0.1, 816], [0.1, 886], [0.1, 953],
     [0.1, 1335], [0.1, 1399], [0.1, 1561], [0.1, 1979], [0.1, 2270],
     [0.1, 2493], [0.1, 2648], [0.1, 2762], [0.1, 2782], [0.1, 2809],
     [0.1, 2846], [0.1, 2910], [0.1, 2960], [0.12, 227], [0.12, 295],
     [0.12, 403], [0.12, 475], [0.12, 537], [0.12, 599], [0.12, 665],
     [0.12, 716], [0.12, 783], [0.12, 835], [0.12, 900], [0.12, 1273],
     [0.12, 1298], [0.12, 1484], [0.12, 1905], [0.12, 2204], [0.12, 2441],
     [0.12, 2606], [0.12, 2733], [0.12, 2753], [0.12, 2783], [0.12, 2824],
     [0.12, 2894], [0.12, 2949], [0.5, 141], [0.5, 203], [0.5, 272],
     [0.5, 322], [0.5, 370], [0.5, 409], [0.5, 449], [0.5, 488], [0.5, 525],
     [0.5, 564], [0.5, 599], [0.5, 759], [0.5, 791], [0.5, 936], [0.5, 1264],
     [0.5, 1583], [0.5, 1875], [0.5, 2116], [0.5, 2329], [0.5, 2369],
     [0.5, 2426], [0.5, 2518], [0.5, 2680], [0.5, 2804], [1.0, 111],
     [1.0, 160], [1.0, 226], [1.0, 270], [1.0, 304], [1.0, 342], [1.0, 377],
     [1.0, 403], [1.0, 438], [1.0, 464], [1.0, 491], [1.0, 578], [1.0, 628],
     [1.0, 766], [1.0, 1009], [1.0, 1271], [1.0, 1543], [1.0, 1801],
     [1.0, 2082], [1.0, 2140], [1.0, 2234], [1.0, 2383], [1.0, 2627],
     [1.0, 2784], [5.0, 129], [5.0, 162], [5.0, 193], [5.0, 220], [5.0, 244],
     [5.0, 265], [5.0, 284], [5.0, 301], [5.0, 322], [5.0, 361], [5.0, 399],
     [5.0, 473], [5.0, 610], [5.0, 760], [5.0, 931], [5.0, 1120], [5.0, 1524],
     [5.0, 1712], [5.0, 2006], [5.0, 2320], [5.0, 2622], [5.0, 2785],
     [10.0, 125], [10.0, 149], [10.0, 172], [10.0, 193], [10.0, 213],
     [10.0, 232], [10.0, 249], [10.0, 265], [10.0, 293], [10.0, 330],
     [10.0, 389], [10.0, 504], [10.0, 634], [10.0, 776], [10.0, 941],
     [10.0, 1289], [10.0, 1556], [10.0, 1997], [10.0, 2322], [10.0, 2624],
     [10.0, 2786]]
)

t_l_points = np.array(
    [[0.001, -5.37], [0.001, -4.715], [0.001, -4.137], [0.001, -3.746],
     [0.001, -3.482], [0.001, -3.273], [0.001, -3.129], [0.001, -2.998],
     [0.001, -2.893], [0.001, -2.811], [0.001, -2.735], [0.001, -2.615],
     [0.001, -2.455], [0.001, -2.28], [0.001, -2.174], [0.001, -1.939],
     [0.001, -1.637], [0.001, -1.604], [0.001, -1.478], [0.001, -1.455],
     [0.001, -1.389], [0.001, -1.373], [0.001, -1.289], [0.001, -1.187],
     [0.005, -6.079], [0.005, -5.522], [0.005, -4.92], [0.005, -4.552],
     [0.005, -4.283], [0.005, -4.061], [0.005, -3.886], [0.005, -3.725],
     [0.005, -3.579], [0.005, -3.46], [0.005, -3.363], [0.005, -3.194],
     [0.005, -2.981], [0.005, -2.661], [0.005, -2.262], [0.005, -2.068],
     [0.005, -1.846], [0.005, -1.902], [0.005, -1.829], [0.005, -1.814],
     [0.005, -1.716], [0.005, -1.76], [0.005, -1.696], [0.005, -1.582],
     [0.01, -6.38], [0.01, -5.87], [0.01, -5.298], [0.01, -4.914],
     [0.01, -4.645], [0.01, -4.429], [0.01, -4.244], [0.01, -4.091],
     [0.01, -3.945], [0.01, -3.82], [0.01, -3.704], [0.01, -3.488],
     [0.01, -3.125], [0.01, -2.709], [0.01, -2.435], [0.01, -2.411],
     [0.01, -2.233], [0.01, -2.188], [0.01, -2.099], [0.01, -2.083],
     [0.01, -2.022], [0.01, -2.016], [0.01, -1.943], [0.01, -1.854],
     [0.05, -7.069], [0.05, -6.587], [0.05, -6.1], [0.05, -5.789],
     [0.05, -5.531], [0.05, -5.336], [0.05, -5.152], [0.05, -4.978],
     [0.05, -4.832], [0.05, -4.706], [0.05, -4.569], [0.05, -4.052],
     [0.05, -3.568], [0.05, -3.642], [0.05, -3.311], [0.05, -3.08],
     [0.05, -2.884], [0.05, -2.755], [0.05, -2.647], [0.05, -2.623],
     [0.05, -2.591], [0.05, -2.551], [0.05, -2.472], [0.05, -2.397],
     [0.1, -7.418], [0.1, -6.957], [0.1, -6.383], [0.1, -6.112], [0.1, -5.88],
     [0.1, -5.686], [0.1, -5.534], [0.1, -5.365], [0.1, -5.246], [0.1, -5.103],
     [0.1, -4.978], [0.1, -4.332], [0.1, -4.281], [0.1, -4.11], [0.1, -3.668],
     [0.1, -3.386], [0.1, -3.167], [0.1, -3.008], [0.1, -2.879], [0.1, -2.856],
     [0.1, -2.821], [0.1, -2.776], [0.1, -2.689], [0.1, -2.617],
     [0.12, -7.526], [0.12, -7.043], [0.12, -6.481], [0.12, -6.183],
     [0.12, -5.97], [0.12, -5.781], [0.12, -5.602], [0.12, -5.477],
     [0.12, -5.326], [0.12, -5.219], [0.12, -5.088], [0.12, -4.432],
     [0.12, -4.436], [0.12, -4.218], [0.12, -3.764], [0.12, -3.472],
     [0.12, -3.244], [0.12, -3.083], [0.12, -2.944], [0.12, -2.922],
     [0.12, -2.886], [0.12, -2.838], [0.12, -2.753], [0.12, -2.673],
     [0.5, -8.415], [0.5, -7.753], [0.5, -7.218], [0.5, -6.913], [0.5, -6.67],
     [0.5, -6.496], [0.5, -6.34], [0.5, -6.2], [0.5, -6.08], [0.5, -5.963],
     [0.5, -5.864], [0.5, -5.447], [0.5, -5.404], [0.5, -5.133], [0.5, -4.636],
     [0.5, -4.255], [0.5, -3.955], [0.5, -3.729], [0.5, -3.534], [0.5, -3.498],
     [0.5, -3.445], [0.5, -3.356], [0.5, -3.189], [0.5, -3.047], [1.0, -8.851],
     [1.0, -8.185], [1.0, -7.56], [1.0, -7.244], [1.0, -7.031], [1.0, -6.831],
     [1.0, -6.664], [1.0, -6.556], [1.0, -6.417], [1.0, -6.325], [1.0, -6.235],
     [1.0, -5.955], [1.0, -5.835], [1.0, -5.514], [1.0, -5.071], [1.0, -4.696],
     [1.0, -4.374], [1.0, -4.106], [1.0, -3.829], [1.0, -3.772], [1.0, -3.679],
     [1.0, -3.527], [1.0, -3.268], [1.0, -3.083], [5.0, -8.57], [5.0, -8.166],
     [5.0, -7.867], [5.0, -7.644], [5.0, -7.469], [5.0, -7.328], [5.0, -7.217],
     [5.0, -7.124], [5.0, -7.015], [5.0, -6.823], [5.0, -6.671], [5.0, -6.401],
     [5.0, -6.008], [5.0, -5.67], [5.0, -5.353], [5.0, -5.058], [5.0, -4.504],
     [5.0, -4.278], [5.0, -3.942], [5.0, -3.603], [5.0, -3.275], [5.0, -3.083],
     [10.0, -8.629], [10.0, -8.325], [10.0, -8.087], [10.0, -7.888],
     [10.0, -7.724], [10.0, -7.584], [10.0, -7.469], [10.0, -7.368],
     [10.0, -7.204], [10.0, -7.016], [10.0, -6.759], [10.0, -6.358],
     [10.0, -6.004], [10.0, -5.695], [10.0, -5.393], [10.0, -4.832],
     [10.0, -4.472], [10.0, -3.954], [10.0, -3.602], [10.0, -3.274],
     [10.0, -3.082]]
)
t_l_points[:, 0] = np.log10(t_l_points[:, 0])

m_values = np.log10(
    [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
     0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.072, 0.075,
     0.08, 0.09, 0.1, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
     0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
     0.072, 0.075, 0.08, 0.09, 0.1, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005,
     0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05,
     0.06, 0.07, 0.072, 0.075, 0.08, 0.09, 0.1, 0.0005, 0.001, 0.002, 0.003,
     0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03,
     0.04, 0.05, 0.06, 0.07, 0.072, 0.075, 0.08, 0.09, 0.1, 0.0005, 0.001,
     0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.012,
     0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.072, 0.075, 0.08, 0.09, 0.1,
     0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
     0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.072, 0.075,
     0.08, 0.09, 0.1, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
     0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
     0.072, 0.075, 0.08, 0.09, 0.1, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005,
     0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05,
     0.06, 0.07, 0.072, 0.075, 0.08, 0.09, 0.1, 0.002, 0.003, 0.004, 0.005,
     0.006, 0.007, 0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05,
     0.06, 0.07, 0.072, 0.075, 0.08, 0.09, 0.1, 0.003, 0.004, 0.005, 0.006,
     0.007, 0.008, 0.009, 0.01, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06,
     0.07, 0.072, 0.075, 0.08, 0.09, 0.1]
)

r_values = np.log10(
    [0.176, 0.166, 0.174, 0.187, 0.2, 0.215, 0.228, 0.242, 0.258, 0.271, 0.285,
     0.307, 0.345, 0.395, 0.408, 0.478, 0.666, 0.665, 0.753, 0.771, 0.847,
     0.84, 0.927, 1.051, 0.148, 0.141, 0.143, 0.148, 0.152, 0.157, 0.161,
     0.166, 0.171, 0.176, 0.181, 0.192, 0.212, 0.261, 0.363, 0.41, 0.51, 0.455,
     0.481, 0.488, 0.55, 0.51, 0.54, 0.611, 0.14, 0.134, 0.136, 0.139, 0.141,
     0.144, 0.147, 0.149, 0.152, 0.155, 0.158, 0.165, 0.192, 0.249, 0.296,
     0.28, 0.323, 0.328, 0.352, 0.357, 0.379, 0.377, 0.401, 0.435, 0.12, 0.122,
     0.124, 0.125, 0.126, 0.126, 0.126, 0.127, 0.127, 0.127, 0.127, 0.137,
     0.153, 0.139, 0.147, 0.158, 0.173, 0.185, 0.198, 0.202, 0.206, 0.212,
     0.225, 0.238, 0.114, 0.117, 0.12, 0.121, 0.122, 0.122, 0.121, 0.121, 0.12,
     0.12, 0.12, 0.129, 0.124, 0.122, 0.126, 0.132, 0.141, 0.15, 0.16, 0.162,
     0.166, 0.17, 0.18, 0.189, 0.113, 0.116, 0.119, 0.12, 0.121, 0.121, 0.12,
     0.12, 0.119, 0.119, 0.119, 0.126, 0.121, 0.119, 0.122, 0.127, 0.135,
     0.142, 0.152, 0.153, 0.157, 0.161, 0.169, 0.178, 0.105, 0.109, 0.112,
     0.113, 0.114, 0.113, 0.113, 0.112, 0.111, 0.11, 0.11, 0.11, 0.107, 0.104,
     0.101, 0.1, 0.101, 0.102, 0.106, 0.107, 0.108, 0.111, 0.119, 0.128, 0.102,
     0.106, 0.109, 0.111, 0.111, 0.11, 0.11, 0.109, 0.108, 0.107, 0.107, 0.106,
     0.103, 0.1, 0.096, 0.093, 0.092, 0.092, 0.094, 0.095, 0.098, 0.102, 0.113,
     0.125, 0.105, 0.105, 0.105, 0.105, 0.104, 0.103, 0.103, 0.102, 0.101, 0.1,
     0.098, 0.095, 0.09, 0.085, 0.082, 0.079, 0.081, 0.083, 0.089, 0.099,
     0.113, 0.125, 0.104, 0.104, 0.103, 0.102, 0.102, 0.101, 0.1, 0.099, 0.098,
     0.096, 0.093, 0.088, 0.083, 0.079, 0.076, 0.078, 0.081, 0.089, 0.099,
     0.113, 0.125]
)

mv_values = np.array(
    [27.57, 25.16, 22.07, 20.51, 19.66, 19.0, 18.48, 17.94, 17.51, 17.17,
     16.85, 16.33, 15.64, 14.88, 14.17, 12.99, 12.21, 11.95, 11.58, 11.51,
     11.44, 11.29, 11.11, 10.92, 29.65, 27.67, 26.0, 24.45, 22.97, 21.84,
     21.06, 20.4, 19.86, 19.43, 19.05, 18.36, 17.38, 15.91, 14.29, 13.26,
     12.56, 12.48, 12.19, 12.14, 11.94, 11.95, 11.73, 11.45, 30.74, 28.5,
     26.99, 25.96, 24.93, 23.89, 22.89, 22.08, 21.33, 20.76, 20.29, 19.46,
     17.94, 16.07, 14.67, 14.16, 13.43, 13.18, 12.83, 12.77, 12.58, 12.53,
     12.26, 11.97, 35.46, 30.95, 28.82, 27.96, 27.4, 26.98, 26.56, 26.12,
     25.68, 25.24, 24.68, 21.96, 19.72, 19.95, 18.27, 16.88, 15.72, 15.0,
     14.46, 14.36, 14.21, 14.02, 13.69, 13.38, 41.98, 32.58, 29.69, 28.71,
     28.09, 27.65, 27.36, 27.03, 26.77, 26.45, 26.1, 23.53, 23.3, 22.3, 19.96,
     18.46, 17.09, 16.08, 15.33, 15.2, 15.01, 14.77, 14.34, 14.02, 43.82,
     33.54, 29.96, 28.93, 28.29, 27.81, 27.48, 27.25, 26.95, 26.71, 26.41,
     24.05, 24.12, 22.95, 20.43, 18.9, 17.51, 16.44, 15.6, 15.47, 15.27, 15.01,
     14.56, 14.19, 56.3, 47.57, 37.05, 32.02, 30.65, 29.6, 29.16, 28.71, 28.4,
     28.14, 27.91, 27.2, 27.11, 26.53, 24.97, 23.11, 21.31, 19.99, 18.84, 18.6,
     18.25, 17.65, 16.54, 15.68, 60.75, 54.15, 44.39, 37.64, 32.62, 31.58,
     30.53, 29.77, 29.37, 29.06, 28.74, 28.09, 27.86, 27.31, 26.4, 25.19,
     23.73, 22.13, 20.44, 20.12, 19.59, 18.67, 16.98, 15.86, 60.05, 55.32,
     50.68, 46.5, 42.71, 39.29, 36.31, 33.73, 33.05, 31.58, 30.2, 29.24, 28.17,
     27.58, 27.09, 26.44, 24.33, 23.11, 21.03, 19.11, 17.02, 15.85, 61.46,
     58.0, 54.62, 51.29, 48.17, 45.19, 42.49, 39.95, 35.23, 32.81, 30.75,
     28.98, 28.18, 27.64, 27.2, 25.69, 24.17, 21.1, 19.1, 17.01, 15.85]
)

t_l_r_interp = LinearNDInterpolator(t_l_points, r_values)
t_l_mv_interp = LinearNDInterpolator(t_l_points, mv_values)
t_l_m_interp = LinearNDInterpolator(t_l_points, m_values)

t_teff_r_interp = LinearNDInterpolator(t_teff_points, r_values)
t_teff_mv_interp = LinearNDInterpolator(t_teff_points, mv_values)
t_teff_m_interp = LinearNDInterpolator(t_teff_points, m_values)
