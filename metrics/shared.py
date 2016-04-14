import numpy as np
import random
import pushover as po
import time
from datetime import timedelta


def set_seeds(seed=1):
    np.random.seed(seed)
    random.seed(seed)


def sec_to_hms(secs):
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    return '{:.02d}:{:02d}:{:05.2f}'.format(h, m, s)


def push_exceptions(fn):
    def ret_fn(*args, **kwargs):
        try:
            start_time = time.time()
            return fn(*args, **kwargs)
        except Exception as e:
            fail_time = time.time()
            exec_time = timedelta(seconds=fail_time - start_time)
            title = "Error in {}:{} after {}".format(fn.__module__, fn.__name__, exec_time)
            message = str(e)

            po.Client().send_message(message, title=title)
            raise e

    return ret_fn


def push_success(fn):
    def ret_fn(*args, **kwargs):
        start_time = time.time()
        result = ret_fn(*args, **kwargs)
        exec_time = timedelta(seconds=time.time() - start_time)
        message = "Completed {}:{} after {}".format(fn.__module__, fn.__name__, exec_time)
        po.Client().send_message(message, title='Success!')
        return result

    return ret_fn
