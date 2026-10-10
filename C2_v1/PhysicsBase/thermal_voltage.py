#热电势
from constants import q,kbv
def thermal_voltage(T):#热电势
    if (T <= 0 ):
        raise ValueError("non_positive Temparature")
    Vt = (kbv * T) / q
    return Vt