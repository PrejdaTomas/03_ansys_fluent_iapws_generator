"""Contains plotting functions for the material properties, as the --plotting argument for the pyIAPWS_toFluent.py is optional, this is also an optional file"""

from typing import Callable
from collections import namedtuple

import matplotlib.pyplot as plt
import matplotlib.figure

plotDescriptor = namedtuple("plotDescriptor", ["VariableName", "xLabel", "yLabel"])
goldenRationer:Callable[[float], tuple[float, float]] = lambda multiplier: (multiplier*1.618, multiplier*1.0)

def genFig2x2(
        title:str,
        size: tuple[float, float] = (4*1.618, 4*1.000)
    ) -> tuple[
        matplotlib.figure.Figure]:
    """Generates the 2 by 2 plot layout using subplots, subplots share the x-axis (e g. Temperature).

    Args:
        - title (str): Title of the figure
        - size (tuple[float, float], optional): x-size and y-size in inches. Defaults to (4*1.618, 4*1.000) - golden ratio.

    Returns:
        matplotlib.figure.Figure: returns the plt.subplots figure
    """
    if not isinstance(title, str):
        raise ValueError("You have not inputted a string input for the plot title: ({0}, {1})".format(
            title,
            title.__class__.__name__
        ))
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(*size)
    fig.suptitle(title)
    return fig
    
def subplotter(
        figure: matplotlib.figure.Figure,
        row: int,
        col: int,
        xData: list[float],
        yData: list[float],
        color: str|tuple[float, float, float],
        variableName: str,
        xLabel: str,
        yLabel: str,
        xCutoff: float = None
    ) -> None:
    """Fills up the subplots of the selected figure

    Args:
        - figure (matplotlib.figure.Figure): selected matplotlib figure
        - row (int): index of the row >= 0 & < 2
        - col (int): index of the column >= 0 & < 2
        - xData (list[float]): 1D iterator with the independent variable values
        - yData (list[float]): 1D iterator with the depedent variable values
        - color (str|tuple[float, float, float]): string representation of the color, but will accept RGB tuple as well
        - variableName (str): name of the plotted dependent variable
        - xLabel (str): x-axis label
        - yLabel (str): y-axis label
        - xCutoff (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see
    """
    axs = figure.axes

    if row == 0 and col == 0: index = 0
    if row == 0 and col == 1: index = 1
    if row == 1 and col == 0: index = 2
    if row == 1 and col == 1: index = 3

    axs[index].plot(xData, yData, color)
    axs[index].set_title(variableName)
    axs[index].set_xlabel(xLabel)
    axs[index].set_ylabel(yLabel)
    if xCutoff is not None:
        axs[index].set_xbound(
                            lower=xData[1] - xCutoff,
                            upper=xData[-2] + xCutoff
                            )    
    axs[index].grid()

def top_left(
    figure: matplotlib.figure.Figure,
    xData: list[float],
    yData: list[float],
    variableName: str,
    xLabel: str,
    yLabel: str,
    xCutoff: float = None
    ) -> None:
    """Fills up the top_left subplot of the selected figure (it is a wrapper for the subplotter function)

    Args:
        - figure (matplotlib.figure.Figure): selected matplotlib figure
        - xData (list[float]): 1D iterator with the independent variable values
        - yData (list[float]): 1D iterator with the depedent variable values
        - variableName (str): name of the plotted dependent variable
        - xLabel (str): x-axis label
        - yLabel (str): y-axis label
        - xCutoff (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see
    """
    subplotter(figure, 0, 0, xData, yData, "black", variableName, xLabel, yLabel, xCutoff)

def top_right(
    figure: matplotlib.figure.Figure,
    xData: list[float],
    yData: list[float],
    variableName: str,
    xLabel: str,
    yLabel: str,
    xCutoff: float = None
    ) -> None:
    """Fills up the top_right subplot of the selected figure (it is a wrapper for the subplotter function)

    Args:
        - figure (matplotlib.figure.Figure): selected matplotlib figure
        - xData (list[float]): 1D iterator with the independent variable values
        - yData (list[float]): 1D iterator with the depedent variable values
        - variableName (str): name of the plotted dependent variable
        - xLabel (str): x-axis label
        - yLabel (str): y-axis label
        - xCutoff (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see
    """
    subplotter(figure, 0, 1, xData, yData, "red", variableName, xLabel, yLabel, xCutoff)
    
def bot_left(
    figure: matplotlib.figure.Figure,
    xData: list[float],
    yData: list[float],
    variableName: str,
    xLabel: str,
    yLabel: str,
    xCutoff: float = None
    ) -> None:
    """Fills up the bot_left subplot of the selected figure (it is a wrapper for the subplotter function)

    Args:
        - figure (matplotlib.figure.Figure): selected matplotlib figure
        - xData (list[float]): 1D iterator with the independent variable values
        - yData (list[float]): 1D iterator with the depedent variable values
        - variableName (str): name of the plotted dependent variable
        - xLabel (str): x-axis label
        - yLabel (str): y-axis label
        - xCutoff (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see
    """
    subplotter(figure, 1, 0, xData, yData, "blue", variableName, xLabel, yLabel, xCutoff)
    
def bot_right(
    figure: matplotlib.figure.Figure,
    xData: list[float],
    yData: list[float],
    variableName: str,
    xLabel: str,
    yLabel: str,
    xCutoff: float
    ) -> None:
    """Fills up the bot_right subplot of the selected figure (it is a wrapper for the subplotter function)

    Args:
        - figure (matplotlib.figure.Figure): selected matplotlib figure
        - xData (list[float]): 1D iterator with the independent variable values
        - yData (list[float]): 1D iterator with the depedent variable values
        - variableName (str): name of the plotted dependent variable
        - xLabel (str): x-axis label
        - yLabel (str): y-axis label
        - xCutoff (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see
    """
    subplotter(figure, 1, 1, xData, yData, "green", variableName, xLabel, yLabel, xCutoff)

def base_layout2x2(
    figTitle: str,
    layoutLabelTups: list[plotDescriptor[str, str, str]],
    independentVariableList: list[float],
    dependentVariableLists: list[list[float]],
    plotSAVEFIG: bool= False,
    plotPATH: str = None,
    plotDPI: int = 250,
    plotSIZE: tuple[float,float] = (4*1.618, 4*1.000),
    plotXCUTOFF: float = None
) -> matplotlib.figure.Figure:
    """
    Args:
        - figTitle (str): title of the figure
        - layoutLabelTups (list[plotDescriptor[str, str, str]]): list of plotDescriptor namedtuples[plotDescriptor("Rho", "T [K]", "Rho [kg/m3]")]
        - independentVariableList (list[float]): x-values
        - dependentVariableLists (list[list[float]]): y-values
        - plotSAVEFIG (bool, optional): if True, saves the figure. Defaults to False.
        - plotPATH (str, optional): path of the saved picture of the plot - do not input the suffix: "C://Users/User/Desktop/figure1". Defaults to None.
        - plotDPI (int, optional): resolution of the picture. Defaults to 250.
        - plotSIZE (tuple[float,float], optional): size of the picture. Defaults to golden ratio (4*1.618, 4*1.000).
        - plotXCUTOFF (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see

    Returns:
        matplotlib.figure.Figure: matplotlib subplot figure containing the data series,
        you do not need to use it, but the function returns the value for further use if needed anyways.
    """
    print("\t-> Generating figure: {0}".format(figTitle))
    fig= genFig2x2(figTitle, plotSIZE)
    top_left(fig, independentVariableList, dependentVariableLists[0], *layoutLabelTups[0], plotXCUTOFF)
    top_right(fig, independentVariableList, dependentVariableLists[1], *layoutLabelTups[1], plotXCUTOFF)
    bot_left(fig, independentVariableList, dependentVariableLists[2], *layoutLabelTups[2], plotXCUTOFF)
    bot_right(fig, independentVariableList, dependentVariableLists[3], *layoutLabelTups[3], plotXCUTOFF)
    if plotSAVEFIG: plt.savefig(plotPATH+ ".png", dpi=plotDPI)
    return fig


def genFig2x1(
        title:str,
        size: tuple[float, float] = (4*1.618, 4*1.000)
    ) -> tuple[
        matplotlib.figure.Figure]:
    """Generates the 2 by 1 plot layout using subplots, subplots share the x-axis (e g. Temperature).

    Args:
        - title (str): Title of the figure
        - size (tuple[float, float], optional): x-size and y-size in inches. Defaults to (4*1.618, 4*1.000) - golden ratio.

    Returns:
        matplotlib.figure.Figure: returns the plt.subplots figure
    """
    if not isinstance(title, str):
        raise ValueError("You have not inputted a string input for the plot title: ({0}, {1})".format(
            title,
            title.__class__.__name__
        ))
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, layout="constrained")
    fig.set_size_inches(*size)
    fig.suptitle(title)
    return fig

def base_layout2x1(
    figTitle: str,
    layoutLabelTups: list[plotDescriptor[str, str, str]],
    independentVariableList: list[float],
    dependentVariableLists: list[list[float]],
    plotSAVEFIG: bool= False,
    plotPATH: str = None,
    plotDPI: int = 250,
    plotSIZE: tuple[float,float] = (4*1.618, 4*1.000),
    plotXCUTOFF: float = None
) -> matplotlib.figure.Figure:
    """
    Args:
        - figTitle (str): title of the figure
        - layoutLabelTups (list[plotDescriptor[str, str, str]]): list of plotDescriptor namedtuples[plotDescriptor("Rho", "T [K]", "Rho [kg/m3]")]
        - independentVariableList (list[float]): x-values
        - dependentVariableLists (list[list[float]]): y-values
        - plotSAVEFIG (bool, optional): if True, saves the figure. Defaults to False.
        - plotPATH (str, optional): path of the saved picture of the plot - do not input the suffix: "C://Users/User/Desktop/figure1". Defaults to None.
        - plotDPI (int, optional): resolution of the picture. Defaults to 250.
        - plotSIZE (tuple[float,float], optional): size of the picture. Defaults to golden ratio (4*1.618, 4*1.000).
        - plotXCUTOFF (float): used for liq and vap values, to which it cuts off the extended value range (to 1 K, to 5000 K), the value is how much of the extension you want to see

    Returns:
        matplotlib.figure.Figure: matplotlib subplot figure containing the data series,
        you do not need to use it, but the function returns the value for further use if needed anyways.
    """
    print("\t-> Generating figure: {0}".format(figTitle))
    fig= genFig2x1(figTitle, plotSIZE)
    top_left(fig,   independentVariableList, dependentVariableLists[0], *layoutLabelTups[0], plotXCUTOFF)
    top_right(fig,  independentVariableList, dependentVariableLists[1], *layoutLabelTups[1], plotXCUTOFF)
    if plotSAVEFIG: plt.savefig(plotPATH+ ".png", dpi=plotDPI)
    return fig
