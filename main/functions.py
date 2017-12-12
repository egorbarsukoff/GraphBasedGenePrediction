import time, datetime


def rev_comp(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def invert(char):
    if char == '+':
        return '-'
    if char == '-':
        return '+'


# noinspection PyDefaultArgument
def log(s, start_time = [0]):
    if s == '__init_log__':
        start_time[0] = time.time()
        return
    print(str(datetime.timedelta(seconds=(time.time()-start_time[0]))) + '  ' + s)