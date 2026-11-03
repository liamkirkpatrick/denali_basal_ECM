"""Public core names and the unchanged source measurement filename prefixes."""


def canonical_core(core):
    return 'DEN13A' if str(core).lower() in ('dic1', 'den13a') else core


def source_core(core):
    return 'dic1' if canonical_core(core) == 'DEN13A' else core


def canonical_stem(stem):
    core, separator, rest = stem.partition('-')
    return canonical_core(core) + separator + rest
