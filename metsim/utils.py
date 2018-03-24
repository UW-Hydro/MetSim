import os
import logging


logger = logging.getLogger(__name__)


def setup_dask(scheduler, num_workers):

    if scheduler is None and num_workers is None:
        return None

    if scheduler is None and num_workers:
        scheduler = 'multiprocessing'  # default when num_workers is not None

    import dask
    options = {}
    if os.path.isfile(scheduler):
        from dask.distributed import Client
        client = Client(scheduler_file=scheduler)
        options['get'] = client.get
        scheduler = 'distributed'
    elif scheduler == 'distributed':
        from dask.distributed import Client
        client = Client(n_workers=num_workers, diagnostics_port=8787)
        logger.info('Dask Client: \n %s \n ' % client)
        options['get'] = client.get
    elif scheduler == 'multiprocessing':
        from multiprocessing import Pool
        import dask.multiprocessing
        options['get'] = dask.multiprocessing.get
        if num_workers is not None:
            options['pool'] = Pool(num_workers)
    elif scheduler == 'threaded':
        from multiprocessing.pool import ThreadPool
        import dask.threaded
        options['get'] = dask.threaded.get
        if num_workers is not None:
            options['pool'] = ThreadPool(num_workers)
    elif scheduler == 'synchronous':
        options['get'] = dask.get
    else:
        raise ValueError('unknown scheduler %s' % scheduler)

    return dask.set_options(**options)
