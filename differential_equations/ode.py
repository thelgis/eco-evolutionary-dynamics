from typing import Tuple
from helpers.type_aliases import DH, DP, Theta, Phi, C, H, P


def predator_prey_no_evolution(time, state: Tuple[H, P], theta: Theta, phi: Phi, c: C) -> Tuple[DH, DP]:
    """
    TODO
    :param time:
    :param state:
    :param theta:
    :param phi:
    :param c:
    :return:
    """

    h, p = state
    dh = h * (1 - theta * h) - (p * h / (1 + h))
    dp = ((phi * p * h) / (1 + h)) - (c * p)

    return dh, dp
