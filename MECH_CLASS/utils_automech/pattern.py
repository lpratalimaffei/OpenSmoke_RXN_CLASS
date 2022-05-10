"""
autoparse.pattern
*****************

Generate re patterns for text parsing.
"""
""" re pattern generators
"""
from re import escape as re_escape



def escape(pattern):
    """ escape special characters in pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return re_escape(pattern)


def maybe(pattern):
    """ a pattern that may or may not be present

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?:{pattern:s})?'


def preceded_by(pattern):
    """ matches if the current position is preceded by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?<={pattern:s})'


def not_preceded_by(pattern):
    """ matches if the current position is not preceded by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?<!{pattern:s})'


def followed_by(pattern):
    """ matches if the current position is followed by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?={pattern:s})'


def not_followed_by(pattern):
    """ matches if the current position is not followed by the pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'(?!{pattern:s})'


def zero_or_more(pattern, greedy=True):
    """ zero or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (rf'(?:{pattern:s})*' if greedy else rf'(?:{pattern:s})*?')


def one_or_more(pattern, greedy=True):
    """ one or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (rf'(?:{pattern:s})+' if greedy else rf'(?:{pattern:s})+?')


def one_of_these(patterns):
    """ any one of a series of patterns

    :param patterns: a series of `re` patterns
    :type patterns: list of strings

    :rtype: str
    """
    patterns_str = '|'.join(patterns)
    return rf'(?:{patterns_str:s})'


def capturing(pattern):
    """ generate a capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return rf'({pattern:s})'


def named_capturing(pattern, name):
    """ generate a named capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param name: a name for the capture
    :type name: str

    :rtype: str
    """
    return rf'(?P<{name:s}>{pattern:s})'


def series(pattern, sep_pattern):
    """ repeated patters with an intervening separator pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param sep_pattern: an `re` pattern
    :type sep_pattern: str

    :rtype: str
    """
    return pattern + zero_or_more(sep_pattern + pattern)


STRING_START = r'\A'
STRING_END = r'\Z'

LINE_START = r'^'
LINE_END = r'$'

WILDCARD = r'[\s\S]'    # literally any character, including spaces
WILDCARD2 = r'.'        # any character except newline
NEWLINE = r'\n'
NONNEWLINE = r'[^\n]'

LINE_FILL = zero_or_more(NONNEWLINE)
LINE = LINE_START + LINE_FILL + LINE_END

SPACE = r'\s'            # space, possibly newline
SPACES = one_or_more(SPACE)
ZSPACES = zero_or_more(SPACE)
LINESPACE = r'[ \t]'     # non-newline space
LINESPACES = one_or_more(LINESPACE)
PADDING = zero_or_more(LINESPACE)
NONSPACE = r'\S'         # any non-space character

UPPERCASE_LETTER = r'[A-Z]'
LOWERCASE_LETTER = r'[a-z]'
LETTER = r'[A-Za-z]'
DIGIT = r'[0-9]'
MINUS = escape('-')
UNDERSCORE = escape('_')
# characters for urlsafe encoding with the `base64` standard library
URLSAFE_CHAR = one_of_these([LETTER, DIGIT, MINUS, UNDERSCORE])

PLUS = escape('+')
SIGN = one_of_these([PLUS, MINUS])
UNSIGNED_INTEGER = one_or_more(DIGIT)
INTEGER = maybe(SIGN) + UNSIGNED_INTEGER

PERIOD = escape('.')
UNSIGNED_FLOAT = one_of_these(
    [zero_or_more(DIGIT) + PERIOD + one_or_more(DIGIT),
     one_or_more(DIGIT) + PERIOD + zero_or_more(DIGIT)])
FLOAT = maybe(SIGN) + UNSIGNED_FLOAT

NUMERICAL_EXPONENT = one_of_these(['E', 'e']) + INTEGER
EXPONENTIAL_INTEGER = INTEGER + NUMERICAL_EXPONENT
EXPONENTIAL_FLOAT = FLOAT + NUMERICAL_EXPONENT

NUMERICAL_EXPONENT_D = one_of_these(['D', 'd']) + INTEGER
EXPONENTIAL_INTEGER_D = INTEGER + NUMERICAL_EXPONENT_D
EXPONENTIAL_FLOAT_D = FLOAT + NUMERICAL_EXPONENT_D

INTEGRAL_NUMBER = one_of_these([EXPONENTIAL_INTEGER, INTEGER])
NUMBER = one_of_these([
    EXPONENTIAL_FLOAT,
    EXPONENTIAL_FLOAT_D,
    EXPONENTIAL_INTEGER,
    FLOAT,
    INTEGER])

UNDERSCORE = escape('_')
VARIABLE_STRING = one_or_more(one_of_these([NONNEWLINE, NONSPACE]))
VARIABLE_NAME = (one_of_these([LETTER, UNDERSCORE]) +
                 one_or_more(one_of_these([LETTER, UNDERSCORE, DIGIT])))


def block_pattern(begin_pattern, end_pattern):
    """ a pattern that will grab all of the block of text
        that exists between the specified begin and end
        patterns
    """
    return (
        begin_pattern +
        capturing(one_or_more(WILDCARD, greedy=False)) +
        end_pattern
    )


def lpadded(pattern, fill_pattern=LINESPACE):
    """ a pattern allowing optional linespaces to the left
    """
    return zero_or_more(fill_pattern) + pattern


def rpadded(pattern, fill_pattern=LINESPACE):
    """ a pattern allowing optional linespaces to the right
    """
    return pattern + zero_or_more(fill_pattern)


def padded(pattern, fill_pattern=LINESPACE):
    """ a pattern allowing optional linespaces to the right
    """
    return zero_or_more(fill_pattern) + pattern + zero_or_more(fill_pattern)


