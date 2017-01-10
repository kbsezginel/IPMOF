# IPMOF Functions for reconstructing unit cell and finding common unit cell
# Date: January 2017
# Author: Kutay B. Sezginel


def lcm(x, tolerance=1):
    """ Calculates least common multiplier for 1 and a given number
        with a tolerance given in percentage """
    tolerance = 100 / tolerance
    multiplier = max(int(x), 1)
    multiplier_found = False
    while(not multiplier_found):
        ratio = round(multiplier / x)
        lower = ratio - ratio / tolerance
        upper = ratio + ratio / tolerance
        if(multiplier / x <= upper and multiplier / x >= lower):
            multiplier_found = True
            break
        else:
            multiplier += 1
    return multiplier
