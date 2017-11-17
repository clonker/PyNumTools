import numpy as _np
import collections as _col


def sliding_window(seq, width=1, fixed_width=False):
    """
    Creates a generator that yields len(seq) tuples of length at most 2*"width"+1 in a sliding window manner along
    the given sequence. If fixed_width=True, the tuples will always have width 2*"width"+1.
    :param seq: the sequence
    :param width: the window width
    :param fixed_width: setting this to True yields tuples that are guaranteed to have length 2*"width"+1
    :return: an iterable
    """
    if len(seq) < 2*width+1:
        raise ValueError("the sequence needs to have at least length 2*width+1 but was {}".format(len(seq)))
    it = iter(seq)
    if not fixed_width:
        win = _col.deque((next(it, None) for _ in range(1 + width)), maxlen=2 * width + 1)
        yield _np.array(win)
        append = win.append
        for e in it:
            append(e)
            yield _np.array(win)
        while len(win) > width + 1:
            win.popleft()
            yield _np.array(win)
    else:
        win = _col.deque((next(it, None) for _ in range(1 + 2*width)), maxlen=2 * width + 1)
        for _ in range(width+1):
            yield _np.array(win)
        append = win.append
        for e in it:
            append(e)
            yield _np.array(win)
        for _ in range(width):
            yield _np.array(win)
