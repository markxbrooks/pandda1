from .base_processors import (
    Processor,
    ProcessorDict,
    ProcessorJoblib,
    ProcessorNotAvailable,
    ProcessWrapper,
)

try:
    from .luigi_processor import ProcessorLuigiSGE
except:
    ProcessorLuigiSGE = ProcessorNotAvailable

basic_processor = Processor()


def wrapper_call(func):
    return func()
