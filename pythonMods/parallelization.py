from twisted.internet import defer, protocol
import os
import subprocess
import pipes
import tempfile


class Scheduler(object):

    def __init__(self, max_threads=float('inf')):
        self.max_threads = max_threads
        self.active_processes = 0
        self.processes = list()
        self.fired_on_completion_deferreds = list()
        self.stop_reactor_deferred = None

    def add_process(self, process):
        self.processes.append(process)
        return process.defer_until_process_completion()

    def run_processes(self):
        from twisted.internet import reactor
        while self.processes and self.active_processes < self.max_threads:
            self._open_next_process()

        reactor.run()

    def _open_next_process(self):
        if self.processes and self.active_processes < self.max_threads:
            p = self.processes.pop(0)
            d = p.launch()
            d.addCallback(self._process_complete_callback)
            self.fired_on_completion_deferreds.append(d)
            self.active_processes += 1
        elif self.stop_reactor_deferred is None and not self.processes:
            # If I get here, the processes are exhausted, and the exit deferred hasn't been created yet
            self.stop_reactor_deferred = defer.DeferredList(self.fired_on_completion_deferreds)
            self.stop_reactor_deferred.addBoth(self._stop_reactor_callback)

    def _stop_reactor_callback(self, results):
        from twisted.internet import reactor
        reactor.stop()

    def _process_complete_callback(self, results):
        self.active_processes -= 1
        self._open_next_process()


class Task(object):

    def __init__(self, action=None):
        self.action = action
        self.fire_on_completion_deferreds = list()

    def launch(self):
        d = self.defer_until_process_completion()
        try:
            result = self.action()
        except:
            # TODO add failure here
            d.errback()
        finally:
            d.callback(result)
        return d

    def fire_fire_on_completion_deferreds(self):
        while self.fire_on_completion_deferreds:
            d = self.fire_on_completion_deferreds.pop()
            try:
                d.callback(None)
            except defer.AlreadyCalledError:
                pass

    def defer_until_process_completion(self):
        d = defer.Deferred()
        self.fire_on_completion_deferreds.append(d)
        return d


class Process(Task, protocol.ProcessProtocol):

    def __init__(self, executable, *args, **kwargs):
        """
        :param executable: Path to the executable to be called
        :param args: Arguments to the Executable being called
        :param output_callback: A function to be called with data from stdout
        :param error_callback: A function to be called with data from stderr
        :param kwargs: Additional keyword arguments, these are largely unused
        """
        super(Process, self).__init__(**kwargs)
        self.exe = executable
        self.args = args
        self.output_callback = kwargs.get('output_callback', echo)
        self.error_callback = kwargs.get('error_callback', echo)
        self.fire_on_completion_deferreds = list()

    def launch(self):
        from twisted.internet import reactor
        reactor.spawnProcess(self, executable=self.exe, args=self.args, env=os.environ)
        return self.defer_until_process_completion()

    def outReceived(self, data):
        self.output_callback(data)

    def errReceived(self, data):
        self.error_callback(data)

    def processEnded(self, reason):
        self.fire_fire_on_completion_deferreds()


class ProcessChain(Task):

    def __init__(self, iterable=None, **kwargs):
        super(ProcessChain, self).__init__(**kwargs)
        if iterable is not None:
            self.processes = [process for process in iterable]
        else:
            self.processes = list()

    def launch(self):
        self._launch_next_process()
        return self.defer_until_process_completion()

    def add_process(self, process):
        self.processes.append(process)
        return process.defer_until_process_completion()

    def get_subprocess_complete_deferreds(self):
        return [process.defer_until_process_completion() for process in self.processes]

    def _launch_next_process(self):
        if self.processes:
            p = self.processes.pop(0)
            d = p.launch()
            d.addCallback(self._launch_next_process_callback)
        else:
            # If I get here, the process chain has been exhausted
            self.fire_fire_on_completion_deferreds()

    def _launch_next_process_callback(self, results):
        self._launch_next_process()


class BatchProcess(Process):
    # TODO: Add memory and time settings

    def __init__(self, *args):
        bsub_path = 'bsub'
        self.submission_args = ['bsub', '-I', '-q', 'interactive']
        self.job_args = list(args)
        Process.__init__(self, bsub_path, *(self.submission_args + self.job_args))

    def launch(self):
        self.args = self.submission_args + self.job_args
        d = Process.launch(self)
        return d

    # TODO: Why doesn't this work on orchestra?
    @staticmethod
    def get_bsub_path():
        proc = subprocess.Popen(['which', 'bsub'], stdout = subprocess.PIPE)
        # TODO: raise process failed error if data(1) is not None
        return proc.communicate()[0]


def echo(stuff):
    print stuff


def process_generator():
    exe = '/bin/echo'
    num = 0
    while num < 15:
        num += 1
        args = ['echo', 'This is process number {}!'.format(num)]
        yield Process(exe, *args)


def test_chains():
    scheduler = Scheduler(max_threads=5)

    print 'Running Chain Processes....'
    for chain in (ProcessChain((process for process in process_generator())) for n in xrange(0, 10)):
        scheduler.add_process(chain)
    scheduler.run_processes()


def test_processes():
    scheduler = Scheduler(max_threads=5)

    print 'Running single processes....'
    for process in process_generator():
        scheduler.add_process(process)
    scheduler.run_processes()


def test_bsub():
    scheduler = Scheduler()

    print 'Submitting Bjobs.....'
    args = ['echo', 'bsub_test', '>', '~/testfile.txt']
    process = BatchProcess(*args)

    scheduler.add_process(process)
    scheduler.run_processes()


def main():
    test_bsub()

if __name__ == '__main__':
    main()
