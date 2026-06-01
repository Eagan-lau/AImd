from extensions.family_runners.ugt_real import ugt_real_runner
from extensions.family_runners.act_real import act_real_runner
from extensions.family_runners.cyp450_real import cyp450_real_runner
from extensions.family_runners.fe2og_real import fe2og_real_runner
from extensions.family_runners.ach_real import ach_real_runner

FAMILY_RUNNERS = {
    "ugt": ugt_real_runner,
    "act": act_real_runner,
    "cyp450": cyp450_real_runner,
    "fe2og": fe2og_real_runner,
    "ach": ach_real_runner,
}
