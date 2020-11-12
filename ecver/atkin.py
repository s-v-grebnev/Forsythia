"""
This module implements interface to call the external ECPP program as threads.
Actually this is not used -- deprecated by more efficient atkin_pro.
"""

import threading
import subprocess as sp

class SingleAtkin(threading.Thread):
    def __init__(self, p, AtkinPath):
        threading.Thread.__init__(self)
        self.p_res = -1
        self.p = p
        self.AtkinPath = AtkinPath

    def run(self):
        with sp.Popen([self.AtkinPath], stdin=sp.PIPE) as fp:
        #                fp.communicate(str(int(self.ui.lineEdit.text(), base = 16))+'\n')
            text = self.p
            text = format('%s' % text)
            text = str(int(text, base=16)) + '\n'
            fp.communicate(bytes(text, 'UTF-8'))
            self.p_res = fp.returncode
#        return self.p_res


def AtkinTest(p, q, AtkinPath):
    thread1 = SingleAtkin(p, AtkinPath)
    thread2 = SingleAtkin(q, AtkinPath)
    thread1.start()
    thread2.start()
    thread1.join()
    thread2.join()
    return (thread1.p_res, thread2.p_res)

