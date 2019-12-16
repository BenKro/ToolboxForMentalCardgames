# ! /usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================================
# gameset.py
#
# 27.07.2019, Benjamin Kromer <benjamin.kromer1@web.de>
#
# @desc: Game information and cards representation, used for playing
#        card games.
# ===================================================================


class GameSet:
    """Game information and cards

    Attributes:
        cards_no (int): number of total game cards
        cards_raw (List[ShortPoints]): integer representing each card forced
            to curve by point multiplication with generator
        cards_masked (List[List[ShortPoints]]): cards_raw masked with pk
        cards_shuffled (List[List[ShortPoints]]): by all players remasked and
            shuffled card deck
        cards_owner (List[str]): information about who owns which card
        cards_unmasked (List[int]): information which element of cards_shuffled
            belongs to which element of cards_raw
        own_id (str): players own ID
        player_no (int): number of players
        player (List[str]): list of player IDs
        player_order (List[int]): order of players for game and shuffle
        player_own_order (int): own order in player_order
        order_share_secret (int): share for generating a player order
        order_share_public (ShortPoint): order_share_secret forced to curve
        order_seed (int): seed calculated by all player for generating the
            player order
    """

    def __init__(self, no_player, no_cards, player, own_id):
        """
        Args:
            no_player (int): number of players in game
            no_cards (int): number of total game cards
            player (List[str]): player IDs
            own_id (str): own ID
        """
        self.cards_no = no_cards

        self.cards_raw = None
        self.cards_masked = None

        self.cards_shuffled = None
        self.cards_owner = [None] * no_cards
        self.cards_unmasked = [None] * no_cards

        self.own_id = own_id
        self.player_no = no_player
        self.player = player
        self.player_order = None
        self.player_own_order = None

        self.order_share_secret = None
        self.order_share_public = None
        self.order_seed = None
