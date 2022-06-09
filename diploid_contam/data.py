from __future__ import annotations

from diploid_contam.models import Interval

EXONS: list[Interval] = [  # hg19
    # GBA
    Interval("1", 155204239, 155204891),
    Interval("1", 155204986, 155205102),
    Interval("1", 155205472, 155205635),
    Interval("1", 155206036, 155206260),
    Interval("1", 155207132, 155207369),
    Interval("1", 155207925, 155208097),
    Interval("1", 155208308, 155208441),
    Interval("1", 155209407, 155209553),
    Interval("1", 155209677, 155209868),
    Interval("1", 155210421, 155210508),
    Interval("1", 155210877, 155211069),
    Interval("1", 155214297, 155214653),
    # CYP21A2
    Interval("6", 32006093, 32006401),
    Interval("6", 32006499, 32006588),
    Interval("6", 32006871, 32007025),
    Interval("6", 32007133, 32007234),
    Interval("6", 32007323, 32007424),
    Interval("6", 32007526, 32007612),
    Interval("6", 32007782, 32007982),
    Interval("6", 32008183, 32008361),
    Interval("6", 32008445, 32008548),
    Interval("6", 32008646, 32009447),
    # PMS2
    Interval("7", 6010556, 6013173),
    Interval("7", 6017219, 6017421),
    Interval("7", 6018227, 6018327),
    Interval("7", 6022455, 6022622),
    Interval("7", 6026390, 6027251),
    Interval("7", 6029431, 6029586),
    Interval("7", 6031604, 6031688),
    Interval("7", 6035165, 6035264),
    Interval("7", 6036957, 6037054),
    Interval("7", 6038739, 6038906),
    Interval("7", 6042084, 6042442),
    Interval("7", 6043321, 6043423),
    Interval("7", 6043603, 6043689),
    Interval("7", 6045523, 6045662),
    Interval("7", 6048438, 6048737),
]
