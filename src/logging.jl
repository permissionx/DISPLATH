#日志记录工具代码
# ===== Professional Logging Utilities =====
using Dates  # 模块导入：导入Julia标准库中的Dates模块，提供时间日期处理功能，用于时间戳生成

# ANSI color codes
const GREEN = "\033[32m"  # 常量定义：绿色文本的ANSI转义序列，转义序列：\033[是ESC字符，32m表示绿色前景色，const关键字：确保常量不可重新赋值
const YELLOW = "\033[33m"  # 黄色文本：ANSI代码33表示黄色，通常用于警告信息
const RED = "\033[31m"  # 红色文本：ANSI代码31表示红色，通常用于错误信息，立即吸引注意力
const BLUE = "\033[34m"  # 蓝色文本：ANSI代码34表示蓝色
const CYAN = "\033[36m"  # 青色文本：ANSI代码36表示青色
const MAGENTA = "\033[35m"  # 洋红色文本：ANSI代码35表示洋红色
const GRAY = "\033[90m"  # 灰色文本：ANSI代码90表示亮黑色（通常显示为灰色），用于次要信息
const WHITE = "\033[37m"  # 白色文本：ANSI代码37表示白色
const RESET = "\033[0m"  # 重置样式：ANSI代码0重置所有文本属性到默认值，必须在使用颜色后重置，避免影响后续输出
const BOLD = "\033[1m"  # 粗体文本：ANSI代码1启用粗体样式，用于强调重要内容
const DIM = "\033[2m"  # 暗淡文本：ANSI代码2启用暗淡（降低亮度）样式，用于次要或调试信息

#文档字符串：Julia的docstring，用于函数文档，参数说明：args...表示可变参数，接受任意数量的参数
"""
    log_info(args...)

Professional logging with timestamp and standard formatting.
"""

function log_info(args...)  # 函数定义：信息级别日志函数，可变参数：使用省略号语法接受任意数量参数
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：获取当前时间并格式化为时:分:秒.毫秒格式，now()来自Dates模块，Dates.format()用于格式化时间字符串
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：使用暗淡的灰色样式输出时间戳，print()不换行输出，RESET重置颜色避免影响后续输出
    
    # Convert all arguments to strings and join them
    message = join(string.(args), "")  # 消息构建：广播string函数到所有参数进行元素级转换，然后使用join连接成单个字符串，广播操作：.符号表示对数组的每个元素应用函数
    print(message)  # 消息输出：输出消息内容，不换行
    println()  # 换行输出：输出换行符，确保每条日志单独一行
end

"""
    log_success(args...)

Success logging with timestamp and green checkmark.
"""
function log_success(args...)  # 函数定义：成功级别日志函数，用于表示操作成功完成
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：与log_info相同，保持时间格式一致性
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：使用标准时间戳格式
    print(GREEN, "✓ ", RESET)  # 成功标志：使用绿色对勾符号，Unicode对勾符号加空格，GREEN设置绿色文本，RESET重置颜色
    
    message = join(string.(args), "")  # 消息构建：与log_info相同的消息构建逻辑
    print(message)  # 消息输出：输出成功消息内容
    println()  # 换行输出：确保日志格式整齐
end

"""
    log_warning(args...)

Warning logging with timestamp and yellow warning symbol.
"""
function log_warning(args...)  # 函数定义：警告级别日志函数，用于需要用户注意的非错误情况
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：标准时间戳格式
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：标准时间戳样式
    print(YELLOW, "⚠ ", RESET)  # 警告标志：使用黄色警告符号，Unicode警告符号加空格，YELLOW设置黄色文本
    
    message = join(string.(args), "")  # 消息构建：将参数转换为字符串并连接
    print(YELLOW, message, RESET)  # 黄色消息：整个消息内容使用黄色显示，强调警告信息的特殊性
    println()  # 换行输出：完成警告日志输出
end

"""
    log_error(args...)

Error logging with timestamp and red X.
"""
function log_error(args...)  # 函数定义：错误级别日志函数，用于表示操作失败或严重问题
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：标准时间戳格式
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：标准时间戳样式
    print(RED, "✗ ", RESET)  # 错误标志：使用红色叉号符号，Unicode叉号符号加空格，RED设置红色文本
    
    message = join(string.(args), "")  # 消息构建：将参数转换为字符串并连接
    print(RED, message, RESET)  # 红色消息：整个错误消息使用红色显示，视觉强调：立即吸引注意力到错误信息
    println()  # 换行输出：完成错误日志输出
end

"""
    log_debug(args...)

Debug logging with timestamp and dim gray color.
"""
function log_debug(args...)  # 函数定义：调试级别日志函数，用于开发调试信息，通常在生产环境中禁用
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：标准时间戳格式
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：标准时间戳样式
    print(DIM, "• ", RESET)  # 调试标志：使用暗淡样式的圆点符号，Unicode圆点符号加空格，DIM设置暗淡样式
    
    message = join(string.(args), "")  # 消息构建：将参数转换为字符串并连接
    print(DIM, message, RESET)  # 暗淡消息：调试信息使用暗淡样式，减少视觉干扰，避免影响主要信息阅读
    println()  # 换行输出：完成调试日志输出
end

"""
    log_section(title)

Creates a section header with consistent formatting.
"""
function log_section(title::String)  # 函数定义：节标题函数，类型注解：明确要求参数为String类型，设计选择：节标题应该是单个连贯的字符串
    timestamp = Dates.format(now(), "HH:MM:SS.sss")  # 时间戳生成：标准时间戳格式
    print(DIM, GRAY, "[", timestamp, "] ", RESET)  # 时间戳输出：标准时间戳样式
    println()  # 换行：在时间戳后换行，使节标题单独成行，增强视觉分离效果
    println(BOLD, "═══ ", title, " ═══", RESET)  # 节标题格式：使用粗体双线装饰，BOLD设置粗体样式，"═══ "和" ═══"为Unicode双线装饰，视觉设计：创建明显的视觉分隔
end

"""
    log_separator()

Creates a simple separator line.
"""
function log_separator()  # 函数定义：分隔线函数，无参数：不需要任何参数，用于创建视觉分隔线
    println(DIM, "─" ^ 60, RESET)  # 分隔线生成：使用暗淡样式的水平线重复60次，"─" ^ 60将Unicode水平线字符重复60次，字符串重复：^操作符在Julia中用于字符串重复
end

# Export the logging functions
export log_info, log_success, log_warning, log_error, log_debug, log_section, log_separator  # 导出语句：使这些函数在模块被导入时可用，模块设计：清晰的API边界，只导出必要的函数