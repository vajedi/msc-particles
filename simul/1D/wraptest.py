def wrap(value, min, max):
    diff = max - min
    value = value.real
    if diff < 0:
        return value
    elif diff == 0:
        return min
    while value > max:
        value -= diff
    while value < min:
        value += diff
    return value
    
    
f = 0.013645462
x = 0.8859728174

f = -0.107708901899
x = 0.89140612528

print 'f+x=', f+x
print wrap(f+x, 0, 1)