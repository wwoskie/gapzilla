class IntervaledFeature:
    def __init__(self, interval, feature):
        self.interval = interval
        self.feature_list = [feature]
        self.feature_num = len(self.feature_list)
        self.feature_lengths = [feature.length]

    def __repr__(self):
        return f"IntervaledFeature(interval={self.interval}, feature_list={self.feature_list}, feature_num={self.feature_num}, feature_lengths={self.feature_lengths})"


class IntervaledGap:
    def __init__(self, interval, features_left, features_right):
        self.interval = interval
        self.features_left = features_left
        self.features_right = features_right

    def __repr__(self):
        return f"IntervaledGap(interval={self.interval}, features_left={self.features_left}, features_right={self.features_right})"


class InsertionSite:
    def __init__(self, start, end, score):
        self.interval = [start, end]
        self.score = score

    def __repr__(self):
        return f"InsertionSite(interval={self.interval}, score={self.score})"
