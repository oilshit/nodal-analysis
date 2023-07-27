from typing import Union

numeric = Union[int, float]

def fahrenheit_to_rankine(T: numeric) -> numeric:
    """
    Converts temperature Fahrenheit into Rankine degree

    input:
        T: int | float
    
    output: int | float
    """
    
    return T + 460

def inch_to_feet(x: numeric) -> numeric:
    return x / 12