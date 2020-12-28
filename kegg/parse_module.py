import re

def fix_module(module):    
    # convert part of a kegg module definition to logical statement for evaluation
    # ----------------------------------------------------------------------------
    # maps [space]:AND, [comma]:OR, [plus]:AND, [minus]:OR
    # adds x[] syntax for evaluation using "f" (see below)
    # example: fix_module('((K00134,K00150) K00927,K11389)') = ((f["K00134"] | f["K00150"]) & f["K00927"] | f["K11389"])
    u = re.sub('\+', ' ', module)    
    u = re.sub('\-', ',', u)
    u = re.sub('([\(, ])([A-Z0-9_])', '\\1f("\\2', u)
    u = re.sub('([A-Z0-9_])([\), ])', '\\1")\\2', u)
    u = re.sub(' ', ' & ', u)
    u = re.sub(',', ' | ', u)
    return(u)

def sum_module(module, x):
    # split kegg module definition into separate components and score them using dictionary
    # -------------------------------------------------------------------------------------
    # module = kegg module definition
    # x = dict that maps each module component to a score (e.g. presence/absence)
    # this function assumes that the module definition is NOT wrapped in parentheses
    # returns a list with the score for each module component
    # example: sum_module(module) = [1, 0, 1, 1]
    def f(a):
        if a == '0':
            return 0
        if a == '1':
            return 1
        if a in x:
            return x[a]
        else:
            return 0    
    while True:
        m = re.findall('(\([^\(\)]*?\))', module)
        if len(m) == 0:
            break
        for mi in m:
            ui = mi
            vi = fix_module(ui)
            module = module.replace(ui, str(eval(vi)))
    out = []
    for mi in module.split():
        out.append(eval(fix_module('(%s)' %(mi))))
    return(out)

