"""
This module implements interface to call the external ECPP program via multiprocesses
"""

import multiprocessing as mp
import subprocess as sp

class SingleAtkin(object):
    def __init__(self, p, AtkinPath):
        self.p_res = -1
        self.p = p
        self.AtkinPath = AtkinPath

    def __call__(self):
        with sp.Popen([self.AtkinPath], stdin=sp.PIPE, stdout=sp.PIPE) as fp:
        #                fp.communicate(str(int(self.ui.lineEdit.text(), base = 16))+'\n')
            text = self.p
            text = format('%s' % text)
            text = str(int(text, base=16)) + '\n'
            fp.communicate(bytes(text, 'UTF-8'))
            self.p_res = fp.returncode
#        return self.p_res


def AtkinTest(p, q, AtkinPath):
    th1 = SingleAtkin(p, AtkinPath)
    th2 = SingleAtkin(q, AtkinPath)
    proc1 = mp.Process(target = th1)
    proc2 = mp.Process(target = th2)

    proc1.start()
    proc2.start()
    proc1.join()
    proc2.join()
    return (th1.p_res, th2.p_res)

